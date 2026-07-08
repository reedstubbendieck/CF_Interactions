# Differential Metabolite Production Underlies Disruption of the Cystic Fibrosis Airway Microbiota by Pathogens

Sydney M. Morabbi<sup>1</sup>, Niladri Bhowmik<sup>1</sup>, Shaz Sutherland<sup>2</sup>, Evelyn A. Wylie<sup>1</sup>, Reagan S. Decker<sup>1</sup>, Akram Al Daerwish<sup>1</sup>, Mercedes Pérez Pérez<sup>1</sup>, Erika I. Lutter<sup>1</sup>, Benjamin Philmus<sup>2</sup>, and Reed M. Stubbendieck<sup>1</sup>

<sup>1</sup>Department of Microbiology and Molecular Genetics, Oklahoma State University, Stillwater, OK 74078
<sup>2</sup>Department of Pharmaceutical Sciences, Oregon State University, Corvallis, Oregon 97331

## Introduction

This repository contains the code necessary to replicate the results and figures of our study of interactions in the cystic fibrosis (CF) airway. 

## Datasets

The code in the Rmd document and the following code snippets expect the raw and derived datasets to be in the `./rawData/` and `./derivedData/` directories, respectively. The data is available here: https://doi.org/10.6084/m9.figshare.32939627.

### Raw datasets

Raw amplicon sequencing reads are available from the NCBI Sequence Read Archive under the run accessions listed in `./rawData/16S_rRNA_gene_amplicon_studies/16S_CF_microbiome.csv`, which also holds the sample metadata used throughout the analysis. Cross-streak bioassay measurements, clinical isolate identifications, LC-MS/MS feature tables, and the paired _Pseudomonas_ genome sequences are in `./rawData/cross-streak/`, `./rawData/isolate_identification/`, `./rawData/lcms/`, and `./rawData/genome_sequences/`, respectively.

### Derived datasets

Processed data are written by the code below and read by the Rmd, including mothur taxonomy summaries (`./derivedData/mothur_output/`) and the comparative genomics outputs (`./derivedData/genome_comparison/`).

## Prerequisites

### Software

*   [antiSMASH](https://docs.antismash.secondarymetabolites.org/install/)
*   [Bakta](https://github.com/oschwengers/bakta)
*   [Bowtie2](https://github.com/BenLangmead/bowtie2)
*   [fastANI](https://github.com/ParBLiSS/FastANI)
*   [fastp](https://github.com/OpenGene/fastp)
*   [Mothur](https://mothur.org/)
*   [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
*   [samtools](https://www.htslib.org/)
*   [snippy](https://github.com/tseemann/snippy)
*   [SRA Toolkit](https://github.com/ncbi/sra-tools)

### R Packages

Statistical analysis was performed in R 4.4.2.

*   [circlize](https://jokergoo.github.io/circlize_book/book/)
*   [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
*   [ComplexUpset](https://cran.r-project.org/web/packages/ComplexUpset/index.html)
*   [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
*   [ggtext](https://cran.r-project.org/web/packages/ggtext/index.html)
*   [multcompView](https://cran.r-project.org/web/packages/multcompView/index.html)
*   [parallelDist](https://cran.r-project.org/web/packages/parallelDist/index.html)
*   [patchwork](https://cran.r-project.org/web/packages/patchwork/index.html)
*   [pvclust](https://cran.r-project.org/web/packages/pvclust/index.html)
*   [readxl](https://cran.r-project.org/web/packages/readxl/index.html)
*   [scales](https://cran.r-project.org/web/packages/scales/index.html)
*   [tidyverse](https://www.tidyverse.org/)
*   [vegan](https://cran.r-project.org/web/packages/vegan/index.html)

## Code

All paths below assume that code is being run from the base project directory.

### 16S rRNA gene amplicon sequencing processing

`download_reads.py` downloads the SRA runs for a metadata file, performs read quality control (fastp trimming, Bowtie2 host-read removal), and classifies the reads with mothur using the batch files in `./scripts/mothur/`. The per-sample `.tax.summary` output is written under `./derivedData/mothur_output/`.

```
conda activate CF_profiling
python3 ./scripts/python/download_reads.py \
  -i ./rawData/16S_rRNA_gene_amplicon_studies/16S_CF_microbiome.csv \
  -o ./derivedData/reads/ \
  --mothur_pe ./scripts/mothur/mothur_batch_file_pe.txt \
  --mothur_se ./scripts/mothur/mothur_batch_file_se.txt
conda deactivate
```

Because the dataset spans many samples, split the run list into chunks and process them in parallel by backgrounding each call:

```
for chunk in a b c d; do
  python3 ./scripts/python/download_reads.py \
    -i runs_${chunk}.txt -o ./derivedData/reads_${chunk}/ \
    --mothur_pe ./scripts/mothur/mothur_batch_file_pe.txt \
    --mothur_se ./scripts/mothur/mothur_batch_file_se.txt &
done
wait
```

Adjust the number of concurrent chunks and the mothur thread counts to suit your hardware. Pass `--local-reads /path/to/fastq/` (or set the `LOCAL_READS` environment variable) to use pre-downloaded FASTQ files instead of downloading from the SRA.

### Pseudomonas comparative genomics

`run_comparison.py` runs fastANI, Bakta, snippy, and antiSMASH on the paired CF _Pseudomonas_ isolates and writes all variant and BGC tables plus the figure plotting data to `./derivedData/genome_comparison/`. Thread count, the Bakta database path, and conda environment names are read from the environment.

```
export BAKTA_DB=/path/to/bakta/db
export THREADS=8
python3 ./scripts/python/run_comparison.py
```

### Analyses

All downstream analyses and figures are produced by `CF_Interactions.Rmd`. Knit the document (`rmarkdown::render("CF_Interactions.Rmd")`) to reproduce them.

## YMLs

*   antismash
*   bakta
*   CF_profiling
*   snippy
