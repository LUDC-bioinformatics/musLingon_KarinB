---
title: "Mouse Lingon diet experiment"
author:
   name: "Shuyi Li & Dmytro Kryvokhyzha"
   email: dmytro.kryvokhyzha@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "02 March, 2021"
output:
    html_document:
      keep_md: true
      toc: true
---




## Project

Mice were fed with different diet for 4 days:

- HFD - high fat diet,
- LFD - low fat diet,
- Lingon - high fat diet with Lingon berries.

HFD vs Lingon is the most interesting comparison. LFD is a control.

RNA-seq data were collected from adipose tissue.

Each group has 5 samples.

**Analysis steps**:

1. QC sequencing data
2. Map and count reads
3. Make pair-wise comparisons between three groups (HFD vs LFD/HFD vs Lingon/LFD vs Lingon) to find the differential expressed (DE) genes
4. Perform enrichment analysis (GO & KEGG)
   - Gene set enrichment analysis
   - Over-representation analysis
5. Check interesing gene, GO terms and KEGG pathways
   * Potential interesting genes
       - FABP4
       - Fsp27
       - CIDEA
       - SLC2A4
       - PNPLA2
       - PLIN1
       - CAV1
       - PPARG
       - DGAT2
       - Anxa2
       - Aacs
       - Acacb
       - Acly
       - Elovls
       - Acots
       - Cidec
       - Insig2
       - ApoA4
   * Interesting GO
      - glucose metabolism
      - adipogenesis
      - mitochondria functions
   * Interesting pathways
      - Lipid synthesis
      - Regulation of lipid metabolism
      - Enzymes in fatty acid activation and oxidation
      - Lipolysis
      - Cholesterol metabolism
      - Glucose uptake

## Data

All the data is located on the Indigo server: `/ludc/Active_Projects/Mouse_Adipocite_Lingonberry/ludc/`

## Prerequisites

You need to install [Conda](https://docs.conda.io/en/latest/) and load the pre-configured conda environment. It should also install all the required programs.


```bash
conda env create -f conf/conda.yml
conda activate LingonProj
```

## Analysis

### QC

Preliminary fastq QC results can be found in: `~/results/tables/multiqc/`

QC reports can be found in: `~/results/reports/`

### Map and count reads

Performed with [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html):


```bash
snakemake -s scr/salmon.smk \
   -j 100 \
   -p --use-conda \
   --cluster-config conf/cluster.yml \
   --cluster "condor_qsub -o logs/{rule}.out -e logs/{rule}.err -l procs={cluster.cores},mem={cluster.ram} -m e -V"
```

Results:

* `results/tables/salmon/{sample}/quant.sf.gz` - transcripts counts
* `results/tables/salmon/{sample}/quant.genes.sf.gz` -gene-level counts

### Differential expression

Performed with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html):


```bash
R -e 'rmarkdown::render("scr/DESeq.Rmd", output_dir="results/reports/")'
```

Results:

* `results/reports/DESeq.html` - notebook describing the analysis
* `results/tables/deseq/LingonProj_DESeqres.csv` - differential expression results with TMP (three results in one table).
* `results/figures/` - plots saved as pdf files

### Enrichment analysis

Performed with [clusterprofiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html):


```bash
# GO analysis
R -e 'rmarkdown::render("scr/GO_analysis.Rmd", output_dir="results/reports/")'
# KEGG analysis
R -e 'rmarkdown::render("scr/KEGG_analysis.Rmd", output_dir="results/reports/")'
```

Results:

* `results/reports/GO_analysis.html` - notebook describing the GO analysis
* `results/tables/GO_results/` - GO analysis results table
* `results/reports/KEGG_analysis.html` - notebook desscribing the KEGG analysis
* `results/tables/KEGG_results/` - KEGG analysis results table

### Select interesting gene

**Note**：interesting_gene.txt is created in advance (case insensitive)


```bash
python scr/select_gene.py -i results/tables/deseq/LingonProj_DESeqres.csv -g data/reference/interesting_gene.txt -o results/tables/deseq/LingonProj_interesting_gene.csv
```

Results:

* `results/tables/deseq/LingonProj_interesting_gene.csv` - DESeq results with only interesting gene (listed above in the analysis step) selected.