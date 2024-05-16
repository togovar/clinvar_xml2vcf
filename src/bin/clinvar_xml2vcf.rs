use clinvar::*;
use flate2::bufread::GzDecoder;
use quick_xml::de::Deserializer;
use quick_xml::events::{BytesStart, Event};
use quick_xml::{Reader, Writer};
use serde::de::Deserialize;
use std::ffi::OsStr;
use std::fs::File;
use std::io::ErrorKind::{AlreadyExists, InvalidInput, NotFound};
use std::io::{self, BufRead, BufReader, BufWriter, Error, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::process::{exit, Command};
use std::str::from_utf8;
use structopt::StructOpt;
use strum::{AsRefStr, EnumString, VariantNames};
use tempfile::tempdir;

const VCF_HEADER: &str = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ClinVar Allele ID">
##INFO=<ID=CONDITIONS,Number=1,Type=String,Description="MedGen:<ID1>/<ID2>/...:<Interpretation1>/<Interpretation2>/...:<SubmissionCount>|MedGen:...">
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=4>
##contig=<ID=5>
##contig=<ID=6>
##contig=<ID=7>
##contig=<ID=8>
##contig=<ID=9>
##contig=<ID=10>
##contig=<ID=11>
##contig=<ID=12>
##contig=<ID=13>
##contig=<ID=14>
##contig=<ID=15>
##contig=<ID=16>
##contig=<ID=17>
##contig=<ID=18>
##contig=<ID=19>
##contig=<ID=20>
##contig=<ID=21>
##contig=<ID=22>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"#;

const EXTENSION_DEBUG_OUTPUT: &'static str = "vcf";
const EXTENSION_OUTPUT: &'static str = "vcf.gz";
const EXTENSION_FAI: &'static str = "gz.fai";
const EXTENSION_GZI: &'static str = "gz.gzi";
const FILE_NAME_TEMP_OUTPUT: &'static str = "output.vcf";
const FILE_NAME_TEMP_SORTED: &'static str = "sorted.vcf.gz";
const FILE_NAME_TEMP_NORMALIZED: &'static str = "normalized.vcf.gz";

#[derive(Debug, EnumString, VariantNames, AsRefStr)]
pub enum Assembly {
    GRCh37,
    GRCh38,
}

#[derive(Debug, StructOpt)]
struct Options {
    /// Just output VCF (do not sort and normalize)
    #[structopt(long)]
    debug: bool,

    /// Overwrite existing file
    #[structopt(long)]
    force: bool,

    /// Continue processing even if an error occurs
    #[structopt(long)]
    ignore_error: bool,

    /// Assembly
    #[structopt(long, possible_values(Assembly::VARIANTS))]
    assembly: Assembly,

    /// Reference fasta
    #[structopt(long, parse(from_os_str))]
    reference: PathBuf,

    /// Path to output
    #[structopt(long, short, parse(from_os_str))]
    output: Option<PathBuf>,

    /// Path to input [*.xml | *.xml.gz]
    #[structopt(parse(from_os_str))]
    input: PathBuf,
}

fn main() -> io::Result<()> {
    let options = Options::from_args();

    if !options.input.exists() {
        Err(Error::new(
            NotFound,
            format!("{}", options.input.to_string_lossy()),
        ))?
    }

    if !options.reference.exists() {
        Err(Error::new(
            NotFound,
            format!("{}", options.reference.to_string_lossy()),
        ))?
    }
    let mut fai = options.reference.clone();
    fai.set_extension(EXTENSION_FAI);
    if !fai.exists() {
        Err(Error::new(NotFound, format!("{}", fai.to_string_lossy())))?
    }
    let mut gzi = options.reference.clone();
    gzi.set_extension(EXTENSION_GZI);
    if !gzi.exists() {
        Err(Error::new(NotFound, format!("{}", gzi.to_string_lossy())))?
    }

    let file_name = options.input.file_name().ok_or(Error::new(
        InvalidInput,
        format!("{}", options.input.to_string_lossy()),
    ))?;

    let output = if let Some(mut o) = options.output {
        if o.is_dir() {
            o.push(file_name);
            o.set_extension(if options.debug {
                EXTENSION_DEBUG_OUTPUT
            } else {
                EXTENSION_OUTPUT
            });
        }
        o
    } else {
        let mut o = Path::new(file_name).to_path_buf();
        o.set_extension(if options.debug {
            EXTENSION_DEBUG_OUTPUT
        } else {
            EXTENSION_OUTPUT
        });
        o
    };

    if output.exists() && !options.force {
        Err(Error::new(
            AlreadyExists,
            format!("{}", output.to_string_lossy()),
        ))?
    }

    let temp_dir = tempdir()?;

    let mut reader = reader_from_path(options.input)?;
    {
        let mut writer = if options.debug {
            BufWriter::new(File::create(&output)?)
        } else {
            BufWriter::new(File::create(temp_dir.path().join(FILE_NAME_TEMP_OUTPUT))?)
        };

        output_vcf(
            &mut reader,
            &mut writer,
            options.assembly.as_ref(),
            options.ignore_error,
        )?;
    }

    if !options.debug {
        if let Err(e) = vcf_sort(
            temp_dir.path().join(FILE_NAME_TEMP_OUTPUT),
            temp_dir.path().join(FILE_NAME_TEMP_SORTED),
        ) {
            std::fs::copy(temp_dir.path().join(FILE_NAME_TEMP_OUTPUT), &output)?;
            eprintln!("Error: {}", e);
            eprintln!("Output temp file to: {}", &output.to_string_lossy());
            exit(1)
        };

        if let Err(e) = vcf_normalize(
            temp_dir.path().join(FILE_NAME_TEMP_SORTED),
            temp_dir.path().join(FILE_NAME_TEMP_NORMALIZED),
            options.reference,
        ) {
            std::fs::copy(temp_dir.path().join(FILE_NAME_TEMP_SORTED), &output)?;
            eprintln!("Error: {}", e);
            eprintln!("Output temp file to: {}", &output.to_string_lossy());
            exit(1)
        };

        std::fs::copy(temp_dir.path().join(FILE_NAME_TEMP_NORMALIZED), &output)?;
        vcf_index(&output)?;
    }

    eprintln!("Output to: {}", &output.to_string_lossy());

    temp_dir.close()
}

fn reader_from_path<T: AsRef<Path>>(path: T) -> io::Result<Reader<Box<dyn BufRead>>> {
    let f = File::open(path.as_ref())?;
    let r: Box<dyn BufRead> = match path.as_ref().extension() {
        Some(ext) if ext == "gz" => Box::new(BufReader::new(GzDecoder::new(BufReader::new(f)))),
        _ => Box::new(BufReader::new(f)),
    };

    Ok(Reader::from_reader(r))
}

fn read_record<R: BufRead>(
    reader: &mut Reader<R>,
    start_tag: &BytesStart,
    buf: &mut Vec<u8>,
) -> Result<Vec<u8>, quick_xml::Error> {
    let tag_name = start_tag.name();
    let mut output_buf = Vec::new();
    let mut w = Writer::new(&mut output_buf);

    w.write_event(Event::Start(start_tag.clone()))?;

    let mut depth = 0;
    loop {
        buf.clear();

        let event = reader.read_event_into(buf)?;

        w.write_event(&event)?;

        match event {
            Event::Start(e) if e.name() == tag_name => depth += 1,
            Event::End(e) if e.name() == tag_name => {
                if depth == 0 {
                    return Ok(output_buf);
                }
                depth -= 1;
            }
            Event::Eof => {
                return Err(quick_xml::Error::UnexpectedEof(
                    "Unexpected end of file (EOF) encountered.".to_string(),
                ));
            }
            _ => {}
        }
    }
}

fn handle_variation_archive(bytes: &Vec<u8>) -> Result<VariationArchive, String> {
    let str = from_utf8(&bytes).map_err(|e| format!("{}", e))?;

    let mut deserializer = Deserializer::from_str(str);

    VariationArchive::deserialize(&mut deserializer).map_err(|e| format!("{}", e))
}

fn output_vcf<R: BufRead, W: Write>(
    mut reader: &mut Reader<R>,
    mut writer: &mut W,
    assembly: &str,
    ignore_error: bool,
) -> io::Result<()> {
    writeln!(writer, "{}", VCF_HEADER)?;

    let mut buf = Vec::new();
    let mut junk_buf = Vec::new();
    loop {
        let event = match reader.read_event_into(&mut buf) {
            Ok(e) => e,
            Err(e) => {
                error!("Error at position {}: {}", reader.buffer_position(), e);
                if ignore_error {
                    continue;
                }
                Err(Error::new(ErrorKind::InvalidData, format!("{}", e)))?
            }
        };

        match event {
            Event::Eof => break,
            Event::Start(start_tag) => match start_tag.name().as_ref() {
                b"VariationArchive" => {
                    match read_record(&mut reader, &start_tag, &mut junk_buf) {
                        Ok(bytes) => match handle_variation_archive(&bytes) {
                            Ok(variant) => output_record(&mut writer, &variant, assembly)?,
                            Err(e) => {
                                error!("Error at position {}: {}", reader.buffer_position(), e);
                                if ignore_error {
                                    continue;
                                }
                                Err(Error::new(ErrorKind::InvalidData, format!("{}", e)))?
                            }
                        },
                        Err(e) => {
                            error!("Error at position {}: {}", reader.buffer_position(), e);
                            if ignore_error {
                                continue;
                            }
                            Err(Error::new(ErrorKind::InvalidData, format!("{}", e)))?
                        }
                    };
                }
                _ => {}
            },
            _ => {}
        }

        buf.clear();
    }

    writer.flush()
}

fn output_record<W: Write>(
    writer: &mut W,
    variant: &VariationArchive,
    assembly: &str,
) -> io::Result<()> {
    if let Some(ref record) = variant.classified_record {
        if let Some(allele) = record.simple_allele.as_ref() {
            if let Some(loc) = extract_location(allele, assembly) {
                let conditions = extract_conditions(record);

                if conditions.is_empty() {
                    warn!(
                        "No ClassifiedCondition associated with MedGen: variation_id = {}",
                        variant.variation_id
                    );
                    return Ok(());
                }

                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t.\t.\tALLELEID={};CONDITIONS={}",
                    loc.0,
                    loc.1,
                    allele.variation_id,
                    loc.2.to_uppercase(),
                    loc.3.to_uppercase(),
                    allele.allele_id,
                    conditions
                )?
            }
        } else {
            warn!(
                "SimpleAllele not found: variation_id = {}",
                variant.variation_id
            );
        }
    } else {
        warn!(
            "ClassifiedRecord not found: variation_id = {}",
            variant.variation_id
        );
        return Ok(());
    }

    Ok(())
}

const DB_MEDGEN: &'static str = "MedGen";

fn extract_conditions(record: &ClassifiedRecord) -> String {
    record
        .rcv_list
        .rcv_accession
        .iter()
        .filter_map(|rcv| {
            let medgen = rcv
                .classified_condition_list
                .classified_condition
                .iter()
                .filter_map(|x| {
                    if x.db.as_ref().map(|x| x.as_str()) == Some(DB_MEDGEN) {
                        x.id.clone()
                    } else {
                        None
                    }
                })
                .collect::<Vec<String>>();

            if medgen.len() != 0 {
                rcv.rcv_classifications
                    .germline_classification
                    .as_ref()
                    .map(|x| {
                        (
                            x.description
                                .text
                                .split(&['/', ';'][..])
                                .map(|x| x.trim().replace(" ", "_").to_lowercase())
                                .collect::<Vec<String>>()
                                .join("/"),
                            x.description.submission_count,
                        )
                    })
                    .map(|(interpretations, submissions)| {
                        format!(
                            "{}:{}:{}:{}",
                            DB_MEDGEN,
                            medgen.join("/"),
                            interpretations,
                            submissions
                        )
                    })
            } else {
                None
            }
        })
        .collect::<Vec<String>>()
        .join("|")
}

fn vcf_sort<T: AsRef<OsStr>>(input: T, output: T) -> io::Result<()> {
    let process = Command::new("bcftools")
        .arg("sort")
        .arg("--output-type")
        .arg("z")
        .arg("--output")
        .arg(output.as_ref())
        .arg(input.as_ref())
        .output()?;

    io::stdout().write_all(&process.stdout)?;
    io::stderr().write_all(&process.stderr)?;

    Ok(())
}

fn vcf_normalize<T: AsRef<OsStr>>(input: T, output: T, reference: T) -> io::Result<()> {
    let process = Command::new("bcftools")
        .arg("norm")
        .arg("--no-version")
        .arg("--output-type")
        .arg("z")
        .arg("--output")
        .arg(output.as_ref())
        .arg("--rm-dup")
        .arg("none")
        .arg("--check-ref")
        .arg("x")
        .arg("--fasta-ref")
        .arg(reference.as_ref())
        .arg(input.as_ref())
        .output()?;

    io::stdout().write_all(&process.stdout)?;
    io::stderr().write_all(&process.stderr)?;

    Ok(())
}

fn vcf_index<T: AsRef<OsStr>>(input: T) -> io::Result<()> {
    let process = Command::new("bcftools")
        .arg("index")
        .arg("--force")
        .arg("--tbi")
        .arg(input.as_ref())
        .output()?;

    io::stdout().write_all(&process.stdout)?;
    io::stderr().write_all(&process.stderr)?;

    Ok(())
}
