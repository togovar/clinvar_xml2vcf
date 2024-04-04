use once_cell::sync::Lazy;
use regex::Regex;
use serde::Deserialize;

#[macro_export]
macro_rules! warn {
    () => {
        eprint!("\n")
    };
    ($($arg:tt)*) => {{
        eprint!("[\x1b[33mWARN\x1b[0m] {}\n", format!($($arg)*))
    }};
}

#[macro_export]
macro_rules! error {
    () => {
        eprint!("\n")
    };
    ($($arg:tt)*) => {{
        eprint!("[\x1b[31mERROR\x1b[0m] {}\n", format!($($arg)*))
    }};
}

#[derive(Debug, Deserialize)]
pub struct VariationArchive {
    #[serde(rename = "@VariationID")]
    pub variation_id: u64,
    #[serde(rename = "@Accession")]
    pub accession: String,
    #[serde(rename = "ClassifiedRecord")]
    pub classified_record: Option<ClassifiedRecord>,
    #[serde(rename = "IncludedRecord")]
    pub included_record: Option<IncludedRecord>,
}

#[derive(Debug, Deserialize)]
pub struct ClassifiedRecord {
    #[serde(rename = "SimpleAllele")]
    pub simple_allele: Option<SimpleAllele>,
    #[serde(rename = "Haplotype")]
    pub haplotype: Option<Haplotype>,
    #[serde(rename = "Genotype")]
    pub genotype: Option<Genotype>,
    #[serde(rename = "RCVList")]
    pub rcv_list: RCVList,
}

#[derive(Debug, Deserialize)]
pub struct IncludedRecord {
    #[serde(rename = "SimpleAllele")]
    pub simple_allele: Option<SimpleAllele>,
    #[serde(rename = "Haplotype")]
    pub haplotype: Option<Haplotype>,
}

#[derive(Debug, Deserialize)]
pub struct Haplotype {
    #[serde(rename = "@VariationID")]
    pub variation_id: u64,
    #[serde(rename = "SimpleAllele")]
    pub simple_allele: Vec<SimpleAllele>,
}

#[derive(Debug, Deserialize)]
pub struct Genotype {
    #[serde(rename = "SimpleAllele")]
    pub simple_allele: Option<Vec<SimpleAllele>>,
    #[serde(rename = "SimpleAllele")]
    pub haplotype: Option<Vec<Haplotype>>,
}

#[derive(Debug, Deserialize)]
pub struct SimpleAllele {
    #[serde(rename = "@AlleleID")]
    pub allele_id: u64,
    #[serde(rename = "@VariationID")]
    pub variation_id: u64,
    #[serde(rename = "Location")]
    pub location: Option<Location>,
}

#[derive(Debug, Deserialize)]
pub struct Location {
    #[serde(default, rename = "SequenceLocation")]
    pub sequence_location: Vec<SequenceLocation>,
}

#[derive(Debug, Deserialize)]
pub struct SequenceLocation {
    #[serde(rename = "@Assembly")]
    pub assembly: String,
    #[serde(rename = "@Chr")]
    pub chr: String,
    #[serde(rename = "@positionVCF")]
    pub pos: Option<u64>,
    #[serde(rename = "@referenceAlleleVCF")]
    pub reference: Option<String>,
    #[serde(rename = "@alternateAlleleVCF")]
    pub alternate: Option<String>,
}

#[derive(Debug, Deserialize)]
pub struct RCVList {
    #[serde(default, rename = "RCVAccession")]
    pub rcv_accession: Vec<RCVAccession>,
}

#[derive(Debug, Deserialize)]
pub struct RCVAccession {
    #[serde(rename = "@Title")]
    pub title: Option<String>,
    #[serde(rename = "@Accession")]
    pub accession: String,
    #[serde(rename = "@Version")]
    pub version: i32,
    #[serde(rename = "ClassifiedConditionList")]
    pub classified_condition_list: ClassifiedConditionList,
    #[serde(rename = "RCVClassifications")]
    pub rcv_classifications: RCVClassifications,
}

#[derive(Debug, Deserialize)]
pub struct ClassifiedConditionList {
    #[serde(default, rename = "ClassifiedCondition")]
    pub classified_condition: Vec<ClassifiedCondition>,
}

#[derive(Debug, Deserialize)]
pub struct ClassifiedCondition {
    #[serde(rename = "@DB")]
    pub db: Option<String>,
    #[serde(rename = "@ID")]
    pub id: Option<String>,
    #[serde(rename = "$text")]
    pub text: String,
}

#[derive(Debug, Deserialize)]
pub struct RCVClassifications {
    #[serde(rename = "GermlineClassification")]
    pub germline_classification: Option<GermlineClassification>,
    #[serde(rename = "SomaticClinicalImpact")]
    pub somatic_clinical_impact: Option<SomaticClinicalImpact>,
    #[serde(rename = "OncogenicityClassification")]
    pub oncogenicity_classification: Option<OncogenicityClassification>,
}

#[derive(Debug, Deserialize)]
pub struct GermlineClassification {
    #[serde(rename = "Description")]
    pub description: Description,
}

#[derive(Debug, Deserialize)]
pub struct SomaticClinicalImpact {
    #[serde(rename = "Description")]
    pub description: Description,
}

#[derive(Debug, Deserialize)]
pub struct OncogenicityClassification {
    #[serde(rename = "Description")]
    pub description: Description,
}

#[derive(Debug, Deserialize)]
pub struct Description {
    #[serde(rename = "@SubmissionCount")]
    pub submission_count: i32,
    #[serde(rename = "$text")]
    pub text: String,
}

pub static REGEX_CHROMOSOME: Lazy<Regex> =
    Lazy::new(|| Regex::new(r"\A([1-9]|1[0-9]|2[0-2]|X|Y|MT)\z").unwrap());
pub static REGEX_ALLELE: Lazy<Regex> = Lazy::new(|| Regex::new(r"\A[ACGT]+\z").unwrap());

/// Extract sequence location from `SimpleAllele`
///
/// # Arguments
///
/// * `allele`: `SimpleAllele`
/// * `assembly`: GRCh38 or GRCh37
///
/// returns: Option<(&String, u64, &String, &String)>
///          (CHROM, POS, REF, ALT)
pub fn extract_location<'a>(
    allele: &'a SimpleAllele,
    assembly: &'a str,
) -> Option<(&'a String, u64, &'a String, &'a String)> {
    allele
        .location
        .as_ref()
        .and_then(|x| x.sequence_location.iter().find(|x| x.assembly == assembly))
        .and_then(|x| match (&x.chr, x.pos, &x.reference, &x.alternate) {
            (c, Some(p), Some(r), Some(a)) => {
                let reference = r.to_uppercase();
                let alternate = a.to_uppercase();

                if !REGEX_CHROMOSOME.is_match(c) {
                    warn!(
                        "Skip chromosome {}: variation_id = {}",
                        c, allele.variation_id
                    );
                    return None;
                }
                if !REGEX_ALLELE.is_match(reference.as_str()) {
                    warn!(
                        "Skip non-ACGT reference: {}, variation_id = {}",
                        reference, allele.variation_id
                    );
                    return None;
                }
                if !REGEX_ALLELE.is_match(alternate.as_str()) {
                    warn!(
                        "Skip non-ACGT alternate: {}, variation_id = {}",
                        alternate, allele.variation_id
                    );
                    return None;
                }
                if reference == alternate {
                    warn!(
                        "Skip ref == alt: {} == {}, variation_id = {}",
                        reference, alternate, allele.variation_id
                    );
                    return None;
                }

                Some((c, p, r, a))
            }
            _ => None,
        })
}
