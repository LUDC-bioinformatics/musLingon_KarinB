---
title: "Mouse Lingon diet experiment"
author:
   name: "Dmytro Kryvokhyzha & Shuyi Li"
   email: dmytro.kryvokhyzha@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "01 december, 2020"
output:
    html_document:
      keep_md: true
      toc: true
---



**NOTE!** We took over this project when most of the analyses have been done.
Unfortunately, there was no clear step-by-step explanation of the analyses.
Below, we re-cover all the analysis steps and describe the analysis.

## Project

Mice were fed with different diet for 4 days:

- HFD - high fat diet,
- LFD - low fat diet,
- Lingon - high fat diet with Lingon berries.

HFD vs Lingon is the most interesting comparison. LFD is a control.

RNA-seq data were collected from adipose tissue.

Each group has 5 samples.

**Analysis steps**:

1. QC sequencing data (NBIS).
2. Obtain the count data  (NBIS).
3. Compare HFD vs LFD vs Lingon (NBIS).
4. Perform enrichment analysis (LUDC-BU).
5. Check interesting GO (LUDC-BU):
  - glucose metabolism
  - adipogenesis
  - mitochondria functions

## Data

All the data is located on the Indigo server: `/ludc/Active_Projects/Mouse_Adipocite_Lingonberry/`

### Original data

The structure of the original data we received is shown in `LingonUppmax_files.txt`.

The backup of all the original data is stored in `b2016057/zipped_stuff.tar.gz`.
Its content is listed in `b2016057/zipped_stuff.tar.gz.content.txt`.

### Current data structure

We cleaned and described all the files to recover the analysis steps. 
We tried to keep the original folder structure and changed the location only of some 
files as described below. We also removed redundant or unnecessary files.

- `Clean` - very raw FASTQ files.

- `nbis_Dorota/` - seems to be a test run that was not used in the final results.

- `nbis_Dorota/Data_Dorota/fastq_Dorota/` - duplicates of `LDF` samples from `Clean`.

- `nbis/Data/fastq` - symlinks to `Clean` with renames sample names.

- `nbis/Data/mm10Genome` - RSEM and other reference files.

- `nbis/Data/mm10Genome/mm10RSEM` and `nbis/Data/mm10Genome/mm10RSEMold` - duplicated RSEM reference (checksums are the same). We deleted them.

- `nbis/Data/mm10Genome/old` - old FASTA and GTF reference files [removed].

- `nbis/Docs` - reads distribution plots (no script, they are not the same as `nbis/Analysis/final/3.3_rseqc/`).

- `nbis/Scripts` - a few scripts (some commands can be recovered from here).

- `nbis/ProjectDesc.txt` - brief description of the project.

- `LingonberryReport/` - uploaded the report and files sent by email (Box.com files).

- `nbis/Scripts/slurm/`- cleaned slurm logs that may be useful to recover the analysis commands.

- `nbis/Analysis/` - originally included several folders with different version of analyses. We classified these folders into:

  -- `nbis/Analysis/final/` - files that are likely were used to produce the final results.

  -- `nbis/Analysis/tests-or-failed/` - failed and test runs that were not used to produce the final results [to remove].

### Final results files

These are the most important files used to produce the final results:

- `nbis/Analysis/final/1.1_fastqc` - FASTQ quality control.

- `nbis/Scripts/rsem-calculate-expression_star_Lingon_HFD_LFD.sbatch` - script to run *RSEM + STAR*.

- `nbis/Analysis/final/2.3.rsem_star` - BAM and reads counts files produced by *RSEM + STAR*.

- `nbis/Analysis/final/rsem_star_slurm_log` - slurm logs of *RSEM + STAR* runs.

- `nbis/Analysis/final/3.3_rseqc` - mapping quality control.

- `nbis/Analysis/final/4.3_multiQC/` - summary of all quality control runs.

-  `/nbis/Analysis/final/EBSeq_2.3.rsem_star/` - *EBSeq* differential expression results (symlink to `nbis/Analysis/final/2.3.rsem_star/Matrix/`).

- `enrichment/` - enrichment results and GO annotation performed by LUDC-BU. See below.

## Analysis (NBIS)

These are the analysis steps which we recovered from the available files and scripts.
Except the GO enrichment analysis, which was done by us
because it was not possible to recover it from the data we received.

### QC sequencing data

Performed with [FastQC v0.7.2](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


```bash
module load 

cd nbis/Analysis/1_fastqc
fastqc Clean/Lingon*/*fq.gz -o nbis/Analysis/1_fastqc/
```

NBIS: Quality looks good, hence no filtering will be performed on the reads before processing.

Results: `nbis/Analysis/final/1.1_fastqc/`.

These and other QC results are summarized with multiQC in `nbis/Analysis/final/4.3_multiQC/multiqc_report.html`
  
### Map and count

Reads were mapped and counted in the [RSEM v1.2.29](http://deweylab.github.io/RSEM/) framework
using [STAR v2.5.1b](https://github.com/alexdobin/STAR) as the aligner.

[Mus_musculus.GRCm38.dna_sm.primary_assembly.fa](ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa) was used as a reference. Only the main 21 scaffolds + MT were used.

Prepare the reference:


```bash
rsem-prepare-reference \
  --gtf nbis/Data/mm10Genome/genes.gtf \
  --star nbis/Data/mm10Genome/genome.fa \
  nbis/Data/mm10Genome/mm10.rsemGenome
#STAR  --runThreadN 1  \
#  --runMode genomeGenerate  \
#  --genomeDir nbis/Data/mm10Genome  \
#  --genomeFastaFiles nbis/Data/mm10Genome/genome.fa  \
#  --sjdbGTFfile nbis/Data/mm10Genome/genes.gtf  \
#  --sjdbOverhang 100  \
#  --outFileNamePrefix nbis/Data/mm10Genome/mm10.rsemGenome
rsem-extract-reference-transcripts \
  nbis/Data/mm10Genome/mm10.rsemGenome 0 \
  nbis/Data/mm10Genome/genes.gtf None 0 \
  nbis/Data/mm10Genome/genome.fa
```
                        
Map and count:


```bash
rsem-calculate-expression --star \
  -p 8 \
  --star-gzipped-read-file \
  --paired-end \
  --estimate-rspd \
  --append-names \
  --output-genome-bam \
  $projDIR/Data/fastq/LFD24_1.fastq.gz $projDIR/Data/fastq/LFD24_2.fastq.gz \
  $projDIR//Data/mm10Genome/mm10.rsemGenome \
  $projDIR/Analysis/2.3.rsem_star/LFD24.rsem 
```

Results:

 - `nbis/Analysis/final/2.3.rsem_star/*.bam` - STAR mapped reads.
 - `nbis/Analysis/final/2.3.rsem_star/*.rsem.genes.results` - gene counts.
 - `nbis/Analysis/final/2.3.rsem_star/*.rsem.isoforms.results` - isoform counts.
 
### Mapping Quality Control

Duplication levels were checked with [RSeQC v2.6.4](https://academic.oup.com/bioinformatics/article/28/16/2184/325191):


```bash
for bamfile in nbis/Analysis/2.3.rsem_star/*transcript.bam; 
  do  
      echo $bamfile;
 read_duplication.py -i ${bamfile} -o ${bamfile}.read-duplication;       
  done
```

Results: `nbis/Analysis/final/3.3_rseqc/`

These and other QC results are summarized with multiQC in `nbis/Analysis/final/4.3_multiQC/multiqc_report.html`

### Compare HFD vs LFD vs Lingon

Differential expression analysis was done in the [RSEM v1.2.29](http://deweylab.github.io/RSEM/) framework
using [EBSeq v1.2.0](https://www.bioconductor.org/packages/release/bioc/html/EBSeq.html):


```bash
rsem-generate-data-matrix  HFD*.rsem.genes.results Lingon*.rsem.genes.results LFD*.rsem.genes.results > Matrix/GeneMat_HFD_Lingon_LDF.tab
rsem-run-ebseq Matrix/GeneMat_HFD_Lingon_LDF.tab 5,5,5 Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults
```

We symlinked `nbis/Analysis/final/2.3.rsem_star/Matrix/` to `/nbis/Analysis/final/EBSeq_2.3.rsem_star/` for convenience.

Results files in `/nbis/Analysis/final/EBSeq_2.3.rsem_star/`:

 - `GeneMat_HFD_Lingon_LDF.Ebseqresults.pattern` - pattern table (see below).
 - `GeneMat_HFD_Lingon_LDF.Ebseqresults` - differential expression results table.
 - `GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab` - only significant differential expression results.
 - `GeneMat_HFD_Lingon_LDF.Ebseqresults.condmeans` - mean expression in each condition.
 - `GeneMat_HFD_Lingon_LDF.Ebseqresults.normalized_data_matrix` - data count table.
 - `pattern2_FDR0.05.list.gz`, `pattern2.list.gz`, `pattern3_FDR0.05.list.gz`, `pattern3.list.gz` - list of all and significant genes in *pattern2* and *pattern3*.
 
Columns in `GeneMat_HFD_Lingon_LDF.Ebseqresults`:

- *Pattern1-Pattern5* - posterior probability of each pattern.
- *MAP* - most likely pattern for a given gene.
- *PPDE* - posterior probability of differential expression for a given gene.

In `GeneMat_HFD_Lingon_LDF.Ebseqresults.pattern` and other files:

- C1 = HDF
- C2 = Lingon
- C3 = LFD

Patterns indicate the following:

- *Pattern2* - LFD specific (genes differ from HFD & Lingon).
- *Pattern3* - Lingon specific (genes differ from HFD & LFD).
- *Pattern4* - HFD specific (genes differ from Lingon and LFD).

### Results

We copied the essential results files to `results/`:

- `results/RSEM/*.results` - gene and isoform counts.
- `results/multiQC` - QC summary.
- `results/EBSeq` - differential expression results.

See the description of the files in the analysis steps.

## Analysis (LUDC-BU)

All analysis related files (code, intermediate files etc.) are located in `enrichment/`.

The main results files are located in `results/enrichment/`. See below.

### Enrichment analysis

We performed the enrichment analysis with two programs: clusterProfiler and gProfiler.
They produce slightly different results: clusterProfiler is more conservative, while
gProfiler tests enrichment in Reactome and WikiPathways in addition to GO and KEGG.

- `results/enrichment/enrichment.html` - report describing the analyses.
- `results/enrichment/clusterprofiler/` - clusterProfiler GO and KEGG enrichment results.
- `results/enrichment/gprofiler/` - gProfiler GO, KEGG, Reactome, and WikiPathways enrichment results.

### Interesting GO terms

Genes related to these functions are of main interested:

  - glucose metabolism
  - adipogenesis
  - mitochondria functions

We have not seen enrichment that can be associated to these functions.
So, we annotated all analysed genes with GO and searched for these function 
among them to see how these genes are expressed.

`results/enrichment/GeneMat_HFD_Lingon_LDF.Ebseqresults.GOannotation.csv` - expression results and GO annotation for all genes.

`results/enrichment/GeneMat_HFD_Lingon_LDF.Ebseqresults.GOannotation_FDR0.05.csv` - subset the differentially expressed genes (FDR=0.05) from this dataset.

`results/enrichment/GeneMat_HFD_Lingon_LDF.Ebseqresults.GOannotation_FDR0.05_glucose-adipo-mitoch.csv` - subset the differentially expressed genes (FDR=0.05) that have terms including words `glucose`, `adipo`, `mitoch`

Column names:

- `HDF Lingon LFD` - mean expression in each condition.
- `*.rsem.genes.results` - the median normalized count for each sample.
- `Pattern1`,	`Pattern2`,	`Pattern3`,	`Pattern4`,	`Pattern5` - posterior probability of the pattern (see above).
- `MAP`  - most likely pattern of differential expression.
- `PPDE` - posterior probability that gene is differentially expressed.
- `Go_id`, `Go_term` - GO annotation of a gene with terms separated by `|`.
- `Go_num` - number of GO terms associated with that gene.


  
