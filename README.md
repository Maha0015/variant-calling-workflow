# Variant Calling Pipeline

##  Overview

This project implements a comprehensive variant calling pipeline that processes raw sequencing reads to identify Single Nucleotide Polymorphisms (SNPs) and insertions/deletions (indels) in human genomic data. The pipeline follows best practices for quality control, read alignment, variant calling, and statistical analysis.

### Key Features
- **Quality Control**: Automated assessment of sequencing read quality using FastQC metrics
- **Read Alignment**: Mapping of sequencing reads to the human reference genome (GRCh38/hg38)
- **Variant Calling**: Identification of SNPs and indels with high-confidence filtering
- **Variant Filtering**: Quality-based filtering (QUAL scores, depth, mapping quality)
- **Statistical Analysis**: Comprehensive variant statistics and quality metrics
- **Visualization**: Publication-quality plots for variant distributions

##  Methods & Workflow

### Pipeline Steps

1. **Data Preprocessing**
   - Input: Raw sequencing reads (FASTQ format)
   - Quality assessment with FastQC
   - Read trimming and adapter removal (if needed)

2. **Read Alignment**
   - Aligner: Minimap2 or BWA-MEM
   - Reference: GRCh38/hg38 human genome
   - Output: Sorted BAM files

3. **Variant Calling**
   - Tool: VarScan mpileup2cns
   - Minimum coverage: 10X
   - Minimum variant frequency: 10%
   - Output: VCF (Variant Call Format)

4. **Quality Filtering**
   - Filter by depth (DP > 50)
   - Filter by quality score (QUAL > 30)
   - Tool: BCFtools

5. **Variant Analysis**
   - Count total variants
   - Classify by type (SNP vs indel)
   - Extract variants in specific genomic regions
   - Analyze variant distribution across chromosomes

6. **Statistical Reporting**
   - Transition/transversion ratios
   - Variant quality distributions
   - Depth of coverage metrics
   - Allele frequency analysis

### Technologies Used
- **Programming**: Python 3.8+, Bash scripting
- **Alignment**: Minimap2, BWA-MEM
- **Variant Calling**: VarScan
- **File Processing**: SAMtools, BCFtools
- **Data Analysis**: Pandas, NumPy
- **Visualization**: Matplotlib, Seaborn

##  Key Results

### Variant Discovery Metrics
- **Total variants called**: Genome-wide SNPs and indels
- **Filtering efficiency**: Reduction in low-quality variants
- **Coverage metrics**: Mean depth, coverage uniformity
- **Variant distribution**: Chromosome-specific variant counts
- **Indel ratio**: Proportion of insertions vs deletions

### Example Analyses Performed
- Extraction of variants in specific genomic regions (e.g., chr16:20800000-30800000)
- Identification of specific variants at known positions
- Counting indels (insertions and deletions)
- Quality-based filtering and comparison

### Quality Metrics
- **QUAL scores**: Phred-scaled quality of variant calls
- **Depth (DP)**: Read coverage at variant positions
- **Mapping quality (MQ)**: Quality of read alignment
- **Base qualities**: Per-base quality scores

##  Usage

### Prerequisites
```bash
pip install -r requirements.txt

# External tools (install separately)
# SAMtools: http://www.htslib.org/
# BCFtools: http://samtools.github.io/bcftools/
# VarScan: http://varscan.sourceforge.net/
```

### Running the Pipeline

```bash
# Interactive analysis in Jupyter
jupyter notebook variant-calling-pipeline.ipynb

# Command-line execution
jupyter nbconvert --execute --to html variant-calling-pipeline.ipynb
```

### Input Data Format
- **FASTQ files**: Raw sequencing reads
- **BAM files**: Aligned reads (sorted and indexed)
- **VCF files**: Variant calls
- **Reference genome**: FASTA format (GRCh38/hg38)

### Example Commands

```bash
# Variant calling with VarScan
varscan mpileup2cns input.mpileup \
  --min-coverage 10 \
  --variants 1 \
  --min-var-freq 0.1 \
  --output-vcf > output.vcf

# Filter variants by depth
bcftools view -i 'FORMAT/DP>50' input.vcf -o filtered.vcf

# Count variants
grep -v "^#" output.vcf | wc -l

# Extract indels
grep -v "^#" output.vcf | awk 'length($4) != length($5)' | wc -l
```

##  Project Structure

```
variant-calling-pipeline/
├── variant-calling-pipeline.ipynb    # Main analysis notebook
├── README.md                          # This file
└── requirements.txt                   # Python dependencies
```

##  Data Sources

### Datasets Used

This analysis was performed on human genome sequencing data:
- **Sample ID**: SRR702076
- **Source**: NCBI Sequence Read Archive (SRA)
- **Type**: Whole genome sequencing (WGS)
- **Reference genome**: GRCh38/hg38

**Data Access:**
- SRA accession: [SRR702076](https://www.ncbi.nlm.nih.gov/sra/?term=SRR702076)
- Can be downloaded using: `fastq-dump SRR702076`

### Data Availability

**To reproduce this analysis:**
1. Download data from SRA: https://www.ncbi.nlm.nih.gov/sra/
2. Download reference genome: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
3. Install required bioinformatics tools (SAMtools, BCFtools, VarScan)
4. Follow pipeline steps in the notebook

##  Biological & Clinical Significance

### Applications of Variant Calling

**1. Precision Medicine**
- Identify disease-causing variants
- Pharmacogenomic variants affecting drug response
- Cancer somatic mutation detection
- Rare disease diagnosis

**2. Population Genetics**
- Genetic diversity studies
- Evolutionary analysis
- Ancestry inference
- Population structure

**3. Clinical Diagnostics**
- Germline variant detection
- Carrier screening
- Prenatal genetic testing
- Cancer genomics (somatic variants)

### Variant Interpretation

Variants identified by this pipeline can be:
- **Pathogenic**: Disease-causing mutations
- **Benign**: Normal genetic variation
- **VUS (Variants of Uncertain Significance)**: Require further investigation
- **Pharmacogenomic**: Affect drug metabolism and response

### Quality Control Importance

Proper variant calling requires:
- Sufficient read depth (typically >10X, ideally >30X)
- High mapping quality (avoid mismapped reads)
- Quality score thresholds (QUAL > 20-30)
- Strand bias assessment
- Allele frequency filtering

##  Dependencies

See `requirements.txt` for Python packages:
```
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
seaborn>=0.12.0
jupyter>=1.0.0
```

**External Tools** (not managed by pip):
- SAMtools 1.15+
- BCFtools 1.15+
- VarScan 2.4+
- Minimap2 or BWA-MEM (for alignment)

##  References & Resources

### Variant Calling Best Practices
- GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651
- Variant calling with SAMtools/BCFtools: http://samtools.github.io/bcftools/howtos/variant-calling.html

### Tools Documentation
- **SAMtools**: http://www.htslib.org/doc/samtools.html
- **BCFtools**: http://samtools.github.io/bcftools/bcftools.html
- **VarScan**: http://varscan.sourceforge.net/using-varscan.html

### File Formats
- **VCF specification**: https://samtools.github.io/hts-specs/VCFv4.2.pdf
- **SAM/BAM specification**: https://samtools.github.io/hts-specs/SAMv1.pd

---
