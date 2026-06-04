# TrioConcordance

**Version:** 1.0.0
**Author:** Xuzhen Li

A Comprehensive Tri-way Genotype Concordance Analyzer.

## Overview

**TrioConcordance** is a robust Python-based bioinformatics tool designed to perform rigorous, locus-by-locus, and sample-by-sample concordance analysis across three independent genotyping datasets (e.g., Whole Genome Sequencing, Target Capture, and a Reference panel).

It implements two main analytical modules:

1. **10-Way State Classification**: Classifies allelic states into 10 mutually exclusive biological categories based on identical-by-state (IBS) logic.
2. **Error Transition Profiling**: Extracts specific mismatch patterns (e.g., Target Capture Errors where `WGS == Reference != Capture`) and profiles the transition frequencies.

## Dependencies

* Python 3.6+
* pandas
* matplotlib

You can install the required packages using pip:

```bash
pip install pandas matplotlib

```

## Data Preparation

The software requires three input genotype matrices (tab-separated). You can extract these matrices directly from your VCF files using `bcftools`.

**Note:** All three input matrices must have identical dimensions and structures (i.e., the same loci and samples in the exact same order).

```bash
# Standard command to extract genotype matrix from VCF
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_wgs.vcf > wgs_matrix.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_cap.vcf > cap_matrix.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' input_ref.vcf > ref_matrix.txt

```

## Usage

```bash
python TrioConcordance.py \
    -w wgs_matrix.txt \
    -c cap_matrix.txt \
    -r ref_matrix.txt \
    -o ./Results \
    -p MyAnalysis

```

### Arguments

| Argument | Description | Required | Default |
| --- | --- | --- | --- |
| `-w`, `--wgs` | Path to WGS genotypes matrix | Yes | - |
| `-c`, `--cap` | Path to Capture genotypes matrix | Yes | - |
| `-r`, `--crr` | Path to Reference genotypes matrix | Yes | - |
| `-o`, `--outdir` | Output directory path | No | `Results` |
| `-p`, `--prefix` | Prefix for output files | No | `Trio` |

## Output Files

The pipeline will generate four standardized output files in the specified directory:

1. **`[prefix]_Global_Stats.tsv`**: A global summary of the 10-way classification across all samples and loci.
2. **`[prefix]_Sample_Stats.tsv`**: Detailed 10-way classification statistics broken down by each individual sample.
3. **`[prefix]_Capture_Error_Transitions.tsv`**: Frequencies of specific genotypic error transitions (profiling Capture errors).
4. **`[prefix]_Concordance_Distribution.pdf`**: A publication-ready stacked bar plot visualizing the proportion of loci in each concordance category per sample.
