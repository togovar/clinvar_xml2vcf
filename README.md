# ClinVar XML to VCF

## Usage

```
clinvar_xml2vcf 0.1.0

USAGE:
    clinvar_xml2vcf [FLAGS] <input> --assembly <assembly> --reference <reference>

FLAGS:
        --debug           Just output VCF (do not sort and normalize)
        --force           Overwrite existing file
    -h, --help            Prints help information
        --ignore-error    Continue processing even if an error occurs
    -V, --version         Prints version information

OPTIONS:
        --assembly <assembly>      Assembly [possible values: GRCh37, GRCh38]
        --reference <reference>    Reference fasta

ARGS:
    <input>    Path to input [*.xml | *.xml.gz]
```

### Prepare sequence references

#### GRCh38

```bash
curl -L "https://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" | gunzip | bgzip > Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bcftools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

#### GRCh37

```bash
curl -L "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz" | gunzip | bgzip > Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
bcftools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
```

### Run

```bash
wget "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz"
clinvar_xml2vcf --ignore-error --assembly GRCh38 --reference Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ClinVarVCVRelease_00-latest.xml.gz 2>&1 | tee log.txt
```
