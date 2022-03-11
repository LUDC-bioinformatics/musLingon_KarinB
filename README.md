---
title: "Mouse Lingon diet experiment"
author:
   name: "Shuyi Li & Dmytro Kryvokhyzha"
   email: dmytro.kryvokhyzha@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "11 mars, 2022"
output:
    html_document:
      keep_md: true
      toc: true
---



## PI

Name: [Karin Stenkula](https://www.ludc.lu.se/karin-stenkula-assistant-professor-pi) & Karin Berger

Email: [karin.stenkula@med.lu.se](mailto:karin.stenkula@med.lu.se)

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
5. Check interesting gene, GO terms and KEGG pathways
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

Performed with [WEB-based GEne SeT AnaLysis Toolkit](http://www.webgestalt.org):


```bash
R -e 'rmarkdown::render("scr/WebGestaltR.Rmd", output_dir="results/reports/")'
```

Results:

* `results/reports/Project_HFD_vs_Lingon_FDR_0_1_ORA` - ORA analysis of HFD_vs_Lingon FDR < 0.1.
* `results/reports/Project_HFD_vs_LFD_FDR_0_01_ORA` - ORA analysis of HFD_vs_LFD FDR < 0.01.
* `results/reports/Project_LFD_vs_Lingon_FDR_0_01_ORA` - ORA analysis of LFD_vs_Lingon FDR < 0.01.

### Select interesting gene

**Note**：interesting_gene.txt is created in advance (case insensitive)


```bash
python scr/select_gene.py -i results/tables/deseq/LingonProj_DESeqres.csv -g data/reference/interesting_gene.txt -o results/tables/deseq/LingonProj_interesting_gene.csv
```

Results:

* `results/tables/deseq/LingonProj_interesting_gene.csv` - DESeq results with only interesting gene (listed above in the analysis step) selected.

### Compare with the published results

We compare our results with the paper [Intact glucose uptake despite deteriorating signaling in adipocytes with high-fat feeding](https://jme.bioscientifica.com/view/journals/jme/60/3/JME-17-0195.xml).

We expect the highest correlation of our results with 4 days results in this paper. 


```bash
R -e 'rmarkdown::render("scr/compare_with_published.Rmd", output_dir="results/reports/")'
```

We see the similarity in fold change and day 4 factor in ANOVA, but mean expression 
is also similar to day 2 and sometimes to day 6 depending at what genes we look.
This may be because the food is a little different - 
the high fat diet with Lingon contained 45% fat, compared with my previous short-term HFD study where the diet contained 58% fat.

## Lab verification

We decided to verify in the lab the following genes:
   
   - mitochondrial fission
   - angiogenesis

These genes are selected in the `scr/WebGestaltR.Rmd` notebook and output to:

* `results/tables/candidates_lab_verification/HFD_vs_Lingon_mitochondrial_fission.csv` - 
DESeq results for genes related to mitochondrial fission.
* `results/tables/candidates_lab_verification/HFD_vs_Lingon_mitochondrial_fission.pdf` -
expression of genes related to mitochondrial fission.
* `results/tables/candidates_lab_verification/HFD_vs_Lingon_angiogenesis.csv` - DESeq results for
genes related to angiogenesis.
* `results/tables/candidates_lab_verification/HFD_vs_Lingon_angiogenesis.pdf` - expression of
genes related to angiogenesis.

## Maqui and Lingon overlap

One of the reviewers requested to check the overlap between gene affected by
maqui berry and our results.

I manually extracted the list of genes from Fig.3 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6769892/ and extracted these genes
from our results:


```bash
head -n 1 results/tables/deseq/LingonProj_DESeqres.csv \
   > results/tables/Maqui_lingon_genes_overlap.csv
grep -f data/Maqui_genes.txt results/tables/deseq/LingonProj_DESeqres.csv \
   >> results/tables/Maqui_lingon_genes_overlap.csv
```

Only *Acaca, Ppargc1a, Acly, Fasn* were among the significant results in the comparison HFD_vs_LFD. These genes are up-regulated in LFD.
*Acaca, Ppargc* are significant in the comparison LFD_vs_Lingon and up-regulated in LFD.

There are no significant matches in HFD_vs_Lingon, although *Acaca* has a p-value of 0.136 and it is up-regulated in Lingon.

So, we may say that we see similar effect only for the *Acaca* gene that is up-regulated by both Maqui and Lingon berries. We can also mention that although we do not get significant difference for these genes in the comparison HFD_vs_Lingon, all these genes except *Prdm16, Cpt1b, Acox3, Prdm16* are expressed on average at higher level in Lingon than in HFD. This overlap is significantly non-random:


```r
dat <- data.frame(
  "Maqui" = c(15, 0),
  "Lingon" = c(11, 4),
  row.names = c("up", "down"),
  stringsAsFactors = FALSE
)
dat
fisher.test(dat, alternative = "greater")
```

So, the effect seems to be partially similar between Maqui and Lingon berries.
