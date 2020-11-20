---
title: "Mouse Lingon diet experiment: Enrichment Analysis"
author:
   name: "Shuyi Li"
   email: shuyi.li@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "20 november, 2020"
output:
  html_document:
    keep_md: true
    toc: true
---



## Libraries


```r
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(gprofiler2)
```

## Seperate genes based on different patterns

Input for `run_select_pattern.sh`: RSEM result (GeneMat_HFD_Lingon_LDF.Ebseqresults/GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab)

Output: Gene list with EnsemblID	and GeneSymbol, classified as diffent patterns
(Output folder: Matrix/pattern)


```bash
#seperate genes based on different patterns
scr/run_select_pattern.sh
#full list of expressed genes
cut -f 1 Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults | \
   sed 's/_/\t/;s/"//g' | \
   tail -n +2 \
   > Matrix/full_gene_list2.csv
```

## Load data


```r
pattern2.FDR0.05 <- read.table("Matrix/pattern/pattern2_FDR0.05.csv",
                               header = TRUE)
pattern3.FDR0.05 <- read.table("Matrix/pattern/pattern3_FDR0.05.csv",
                               header = TRUE)
pattern4.FDR0.05 <- read.table("Matrix/pattern/pattern4_FDR0.05.csv",
                               header = TRUE)
pattern_table <- read.table(gzfile("Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults.pattern.gz"),
                            header = TRUE)
full_gene <- read.table("Matrix/full_gene_list.csv", header = FALSE)
names(full_gene) <- c("EnsemblID","GeneSymbol")
```

## RSEM result summary

Summarize the number of genes classified to each pattern (FDR cutoff = 0.05)


```r
names(pattern_table) <- c("HDF","Lingon","LFD")
RSEM_summary <-data.frame(c(0,155,40,123,0,318),
                 row.names = c("pattern1.FDR0.05","pattern2.FDR0.05",
                               "pattern3.FDR0.05","pattern4.FDR0.05",
                               "pattern5.FDR0.05","Sum"))
names(RSEM_summary) <- c("number of genes")
knitr::kable(pattern_table)
```



|         | HDF| Lingon| LFD|
|:--------|---:|------:|---:|
|Pattern1 |   1|      1|   1|
|Pattern2 |   1|      1|   2|
|Pattern3 |   1|      2|   1|
|Pattern4 |   1|      2|   2|
|Pattern5 |   1|      2|   3|

```r
knitr::kable(RSEM_summary)
```



|                 | number of genes|
|:----------------|---------------:|
|pattern1.FDR0.05 |               0|
|pattern2.FDR0.05 |             155|
|pattern3.FDR0.05 |              40|
|pattern4.FDR0.05 |             123|
|pattern5.FDR0.05 |               0|
|Sum              |             318|

Explanation:

- **Pattern1** - no differential change between any of the experiments.
- **Pattern2** - differentially expressed (DE) in LFD compared to HFD & Lingon.
- **Pattern3** - DE in Lingon compared to LFD & HFD.
- **Pattern4** - DE in HFD compared to Lingon & LFD.

## ClusterProfiler

[ClusterProfiler](https://guangchuangyu.github.io/software/clusterProfiler)


```r
ego <- function(pattern_data){
   ego_pattern <- enrichGO(gene = pattern_data[,1], # test gene list
                           universe = full_gene[,1], # background gene list
                           OrgDb = org.Mm.eg.db, 
                           keyType = "ENSEMBL", 
                           ont = 'ALL',       
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
   file_name <- paste("results/clusterprofiler/",
                      deparse(substitute(pattern_data)),
                      "_GO_enrich.csv", sep = "")
   write.table(as.data.frame(ego_pattern), file_name, row.names = F, sep = "\t")
   return(ego_pattern)
}

full_entrez <- bitr(full_gene[,1],
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
eKEGG <- function(pattern_data){
   pattern_entrez <- bitr(pattern_data[,1],
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Mm.eg.db)
   eKEGG_pattern <- enrichKEGG(gene = pattern_entrez[,2], # test gene list
                           universe = full_entrez[,2], # background gene list
                           organism='mmu',
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
   eKEGG_pattern_genes <- setReadable(eKEGG_pattern,
                                      OrgDb = org.Mm.eg.db,
                                      keyType="ENTREZID")
   file_name <- paste("results/clusterprofiler/",
                      deparse(substitute(pattern_data)),
                      "_KEGG_enrich.csv", sep = "")
   write.table(as.data.frame(eKEGG_pattern_genes),
             file_name, row.names = F, sep = "\t")
   return(eKEGG_pattern_genes)
}
```

### GO pattern2.FDR0.05 (LFD specific)


```r
ego_pattern2.FDR0.05 <- ego(pattern2.FDR0.05)
ego_pattern2.FDR0.05_n  <- dim(ego_pattern2.FDR0.05)[1]
ego_pattern2.FDR0.05_n
```

```
## [1] 30
```

Gene counts and significance of enrichment:


```r
barplot(ego_pattern2.FDR0.05, showCategory=ego_pattern2.FDR0.05_n)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Gene ratio, counts and significance of enrichment:


```r
dotplot(ego_pattern2.FDR0.05, showCategory=ego_pattern2.FDR0.05_n)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Differentially expressed genes that belong enriched terms:


```r
heatplot(ego_pattern2.FDR0.05)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

### GO pattern3.FDR0.05 (Lingon specific)


```r
ego_pattern3.FDR0.05 <- ego(pattern3.FDR0.05)
ego_pattern3.FDR0.05_n  <- dim(ego_pattern3.FDR0.05)[1]
ego_pattern3.FDR0.05_n
```

```
## [1] 10
```

Gene counts and significance of enrichment:


```r
barplot(ego_pattern3.FDR0.05,showCategory=ego_pattern3.FDR0.05_n,drop=T)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

Gene ratio, counts and significance of enrichment:


```r
dotplot(ego_pattern3.FDR0.05,showCategory=ego_pattern3.FDR0.05_n)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

Differentially expressed genes that belong enriched terms:


```r
heatplot(ego_pattern3.FDR0.05)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

### GO pattern4.FDR0.05 (HFD specific)


```r
ego_pattern4.FDR0.05 <- ego(pattern4.FDR0.05)
ego_pattern4.FDR0.05_n  <- dim(ego_pattern4.FDR0.05)[1]
ego_pattern4.FDR0.05_n
```

```
## [1] 27
```

Gene counts and significance of enrichment:


```r
barplot(ego_pattern4.FDR0.05,showCategory=ego_pattern4.FDR0.05_n,drop=T)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

Gene ratio, counts and significance of enrichment:


```r
dotplot(ego_pattern4.FDR0.05,showCategory=ego_pattern4.FDR0.05_n)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Differentially expressed genes that belong enriched terms:


```r
heatplot(ego_pattern4.FDR0.05)
```

![](/mnt/Performance/Science/musLingon_KarinB/results/enrichment/enrichment_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


### KEGG pattern2.FDR0.05 (LFD specific)


```r
eKEGG_pattern2.FDR0.05 <- eKEGG(pattern2.FDR0.05)
knitr::kable(as.data.frame(eKEGG_pattern2.FDR0.05))
```



|         |ID       |Description                         |GeneRatio |BgRatio | pvalue| p.adjust|  qvalue|geneID                                                   | Count|
|:--------|:--------|:-----------------------------------|:---------|:-------|------:|--------:|-------:|:--------------------------------------------------------|-----:|
|mmu04610 |mmu04610 |Complement and coagulation cascades |10/68     |88/6987 |      0|  1.6e-06| 1.6e-06|Serpine1/Fga/Kng1/Fgg/Fgb/Serpina1e/Serpina1c/Plg/F2/F2r |    10|

### KEGG pattern3.FDR0.05 (Lingon specific)


```r
eKEGG_pattern3.FDR0.05 <- eKEGG(pattern3.FDR0.05)
knitr::kable(as.data.frame(eKEGG_pattern3.FDR0.05))
```



|ID |Description |GeneRatio |BgRatio | pvalue| p.adjust| qvalue|geneID | Count|
|:--|:-----------|:---------|:-------|------:|--------:|------:|:------|-----:|

### KEGG pattern4.FDR0.05 (HFD specific)


```r
eKEGG_pattern4.FDR0.05 <- eKEGG(pattern4.FDR0.05)
knitr::kable(as.data.frame(eKEGG_pattern4.FDR0.05))
```



|         |ID       |Description            |GeneRatio |BgRatio  |   pvalue|  p.adjust|    qvalue|geneID                                                 | Count|
|:--------|:--------|:----------------------|:---------|:--------|--------:|---------:|---------:|:------------------------------------------------------|-----:|
|mmu04015 |mmu04015 |Rap1 signaling pathway |9/51      |201/6987 | 1.20e-05| 0.0022877| 0.0016604|Prkd2/Mapk3/Calml3/Arap3/Fgfr4/Pik3r3/Adcy8/Cdh1/Gnai2 |     9|
|mmu04024 |mmu04024 |cAMP signaling pathway |8/51      |201/6987 | 8.97e-05| 0.0085247| 0.0061869|Mapk3/Calml3/Arap3/Pde4b/Pik3r3/Ppp1r1b/Adcy8/Gnai2    |     8|

## gProfiler

[gProfiler2](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)


```r
namesP <- c("term_id", "source", "p_value", 'term_size', 'intersection_size',
            "term_name","intersection")
ggo <- function(pattern_data){
   gost_pattern <- gost(pattern_data[,2],
                        domain_scope = "custom",
                        custom_bg = full_gene[,2],
                        sources = c('GO', 'KEGG', 'REAC', 'WP'),
                        organism ='mmusculus',
                        significant = TRUE,
                        evcodes=TRUE)
   file_name <- paste("results/gprofiler/",
                      deparse(substitute(pattern_data)),
                      "_gProfiler_enrich.csv", sep = "")
   write.table(as.data.frame(gost_pattern$result[,namesP]),
            file_name, row.names = F, sep = "\t")
   return(gost_pattern)
}
```

### pattern2.FDR0.05 (LFD specific)


```r
gost_pattern2 <- ggo(pattern2.FDR0.05)
gostplot(gost_pattern2)
```

<!--html_preserve--><div id="htmlwidget-f7a43ec50f94fd85160a" style="width:960px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-f7a43ec50f94fd85160a">{"x":{"data":[{"x":[137.842136997783,101.947433056461,82.7073785951368,124.059353920607,137.845763091726,124.073858296379,125.832513858713,113.315237567629,141.758318456176,146.135013845324,85.4233229584108,125.847018234485,132.087525910312,98.0530081617267,125.915914019401,154.649082423384,93.4551210420584,82.7110046890797,98.0457559738407,154.645456329441,125.912287925458,99.9494552938927,124.055727826664,127.000116108345,112.234661572628,134.83247902513,101.951059150404,87.1565958631438,116.158095218906,79.9225384469465,116.154469124963,75.6256171245436,161.560417478658,85.2021312278905,126.155236219636,164.678858269601,125.172564761095,127.029124859888,114.32329168377,177.185256278856,127.192299087321,105.145647914149,127.221307838865,159.939553486157,81.3874803999008,88.0776237246546,95.3696986439392,159.975814425587,87.1529697692009,116.21248662805,134.817974649358,113.021523958249,85.2057573218335,102.955487172603,127.003742202288,87.4031702512648,100.493369385336,85.5466101524713,160.349302101711,124.255162993527,127.206803463093,85.2093834157764,86.7831081870193,161.556791384716,167.91696016066,82.442673737301,88.0522410670539,82.4535520191298,125.908661831515,174.35690300335,87.2363699298888,86.7867342809623,164.675232175658,85.2130095097194,94.5248187552304,88.0812498185976,124.25878908747,114.127482610851,143.756296218744,162.651871755488,105.138395726263,154.652708517327,98.0493820677837,124.03759735695,104.953464935172,177.181630184913,99.9530813878356,77.6308470749983,139.673314438976,111.618225602325,127.210429557036],"y":[8.84323075353602,7.32077794879062,7.07929550622592,7.07124620263153,6.99294644861301,6.35136907959912,6.05905570227234,5.48383063827338,5.35900512962933,5.26250398246404,4.80773952166713,4.6192236703132,4.61280741435517,4.58963070903354,4.51859669152703,4.51859669152703,4.45026050038976,4.43703583188453,4.04045398699056,3.99711104468477,3.95436335019714,3.76866850213384,3.70524114485024,3.61306965707427,3.41500161339685,3.33430140006676,3.32007942932627,3.22968815630754,3.18265915915966,2.93928001400753,2.82358601432427,2.78908721122829,2.75258354958184,2.73471205011197,2.69649616065186,2.68399533823032,2.65981474966087,2.62293471902447,2.62277865543906,2.5958378601839,2.58187107723025,2.52207359927266,2.50676134412988,2.48359390838005,2.42351864759007,2.4022098088888,2.37732232262891,2.28118257909338,2.28038062762357,2.2004426694673,2.19898888954856,2.19620678019172,2.08033243632169,2.07586069330702,2.00645340171313,1.92875987794647,1.85062536876228,1.83327822290436,1.82795627215623,1.71927115439975,1.71697118839139,1.68528122906427,1.66931380008923,1.64964995050894,1.64964995050894,1.64321151497225,1.63945034026202,1.60472144484074,1.60472144484074,1.58567702724672,1.54643152188132,1.52590051023133,1.45596160540875,1.44705868917888,1.44667209795253,1.44499361798045,1.43744344036463,1.41680379508109,1.39985928678837,1.38564320500443,1.37467471634444,1.37467471634444,1.37467471634444,1.37169200974724,1.37163080795968,1.34491312889214,1.33248589423044,1.33217502268124,1.32674534083869,1.32446405164023,1.3234021132404],"text":["GO:0065007 (10242) <br> biological regulation <br> 1.435e-09","GO:0032501 (6189) <br> multicellular organismal process <br> 4.778e-08","GO:0008150 (16263) <br> biological_process <br> 8.331e-08","GO:0048519 (4885) <br> negative regulation of biological process <br> 8.487e-08","GO:0065008 (3540) <br> regulation of biological quality <br> 1.016e-07","GO:0048523 (4468) <br> negative regulation of cellular process <br> 4.453e-07","GO:0050789 (9684) <br> regulation of biological process <br> 8.729e-07","GO:0042730 (19) <br> fibrinolysis <br> 3.282e-06","GO:0071704 (9613) <br> organic substance metabolic process <br> 4.375e-06","GO:0080090 (5229) <br> regulation of primary metabolic process <br> 5.464e-06","GO:0009987 (14055) <br> cellular process <br> 1.557e-05","GO:0050794 (9266) <br> regulation of cellular process <br> 2.403e-05","GO:0060255 (5540) <br> regulation of macromolecule metabolic process <br> 2.439e-05","GO:0030195 (44) <br> negative regulation of blood coagulation <br> 2.573e-05","GO:0050819 (45) <br> negative regulation of coagulation <br> 3.030e-05","GO:1900047 (45) <br> negative regulation of hemostasis <br> 3.030e-05","GO:0019222 (5973) <br> regulation of metabolic process <br> 3.546e-05","GO:0008152 (10097) <br> metabolic process <br> 3.656e-05","GO:0030193 (79) <br> regulation of blood coagulation <br> 9.111e-05","GO:1900046 (80) <br> regulation of hemostasis <br> 1.007e-04","GO:0050818 (81) <br> regulation of coagulation <br> 1.111e-04","GO:0031323 (5418) <br> regulation of cellular metabolic process <br> 1.703e-04","GO:0048518 (5565) <br> positive regulation of biological process <br> 1.971e-04","GO:0051171 (5090) <br> regulation of nitrogen compound metabolic process <br> 2.437e-04","GO:0042221 (3826) <br> response to chemical <br> 3.846e-04","GO:0061045 (66) <br> negative regulation of wound healing <br> 4.631e-04","GO:0032502 (5562) <br> developmental process <br> 4.785e-04","GO:0010605 (2548) <br> negative regulation of macromolecule metabolic process <br> 5.893e-04","GO:0044238 (9061) <br> primary metabolic process <br> 6.567e-04","GO:0006807 (8571) <br> nitrogen compound metabolic process <br> 1.150e-03","GO:0044237 (9206) <br> cellular metabolic process <br> 1.501e-03","GO:0003008 (1547) <br> system process <br> 1.625e-03","GO:1902042 (28) <br> negative regulation of extrinsic apoptotic signaling pathway via death domain receptors <br> 1.768e-03","GO:0009892 (2772) <br> negative regulation of metabolic process <br> 1.842e-03","GO:0050896 (7388) <br> response to stimulus <br> 2.011e-03","GO:1903035 (82) <br> negative regulation of response to wounding <br> 2.070e-03","GO:0048856 (5147) <br> anatomical structure development <br> 2.189e-03","GO:0051179 (5542) <br> localization <br> 2.383e-03","GO:0043170 (8265) <br> macromolecule metabolic process <br> 2.384e-03","GO:2000352 (30) <br> negative regulation of endothelial cell apoptotic process <br> 2.536e-03","GO:0051235 (320) <br> maintenance of location <br> 2.619e-03","GO:0034116 (14) <br> positive regulation of heterotypic cell-cell adhesion <br> 3.006e-03","GO:0051246 (2613) <br> regulation of protein metabolic process <br> 3.113e-03","GO:1901564 (5851) <br> organonitrogen compound metabolic process <br> 3.284e-03","GO:0007275 (4726) <br> multicellular organism development <br> 3.771e-03","GO:0010883 (58) <br> regulation of lipid storage <br> 3.961e-03","GO:0019915 (91) <br> lipid storage <br> 4.194e-03","GO:1901575 (1838) <br> organic substance catabolic process <br> 5.234e-03","GO:0010604 (3225) <br> positive regulation of macromolecule metabolic process <br> 5.243e-03","GO:0044262 (289) <br> cellular carbohydrate metabolic process <br> 6.303e-03","GO:0061041 (137) <br> regulation of wound healing <br> 6.324e-03","GO:0042592 (1757) <br> homeostatic process <br> 6.365e-03","GO:0009893 (3513) <br> positive regulation of metabolic process <br> 8.311e-03","GO:0032879 (2614) <br> regulation of localization <br> 8.397e-03","GO:0051172 (2205) <br> negative regulation of nitrogen compound metabolic process <br> 9.853e-03","GO:0010675 (149) <br> regulation of cellular carbohydrate metabolic process <br> 1.178e-02","GO:0031639 (20) <br> plasminogen activation <br> 1.411e-02","GO:0010033 (2900) <br> response to organic substance <br> 1.468e-02","GO:1901700 (1637) <br> response to oxygen-containing compound <br> 1.486e-02","GO:0048583 (3507) <br> regulation of response to stimulus <br> 1.909e-02","GO:0051239 (2933) <br> regulation of multicellular organismal process <br> 1.919e-02","GO:0009894 (843) <br> regulation of catabolic process <br> 2.064e-02","GO:0010467 (4991) <br> gene expression <br> 2.141e-02","GO:1902041 (46) <br> regulation of extrinsic apoptotic signaling pathway via death domain receptors <br> 2.241e-02","GO:1904036 (46) <br> negative regulation of epithelial cell apoptotic process <br> 2.241e-02","GO:0007596 (163) <br> blood coagulation <br> 2.274e-02","GO:0010876 (396) <br> lipid localization <br> 2.294e-02","GO:0007599 (165) <br> hemostasis <br> 2.485e-02","GO:0050817 (165) <br> coagulation <br> 2.485e-02","GO:1905952 (166) <br> regulation of lipid localization <br> 2.596e-02","GO:0010628 (2107) <br> positive regulation of gene expression <br> 2.842e-02","GO:0010468 (4048) <br> regulation of gene expression <br> 2.979e-02","GO:1903034 (173) <br> regulation of response to wounding <br> 3.500e-02","GO:0009895 (285) <br> negative regulation of catabolic process <br> 3.572e-02","GO:0019538 (5075) <br> protein metabolic process <br> 3.575e-02","GO:0010884 (25) <br> positive regulation of lipid storage <br> 3.589e-02","GO:0048584 (2029) <br> positive regulation of response to stimulus <br> 3.652e-02","GO:0043086 (717) <br> negative regulation of catalytic activity <br> 3.830e-02","GO:0072378 (9) <br> blood coagulation, fibrin clot formation <br> 3.982e-02","GO:1902373 (52) <br> negative regulation of mRNA catabolic process <br> 4.115e-02","GO:0034114 (26) <br> regulation of heterotypic cell-cell adhesion <br> 4.220e-02","GO:1900048 (26) <br> positive regulation of hemostasis <br> 4.220e-02","GO:0030194 (26) <br> positive regulation of blood coagulation <br> 4.220e-02","GO:0048513 (3150) <br> animal organ development <br> 4.249e-02","GO:0033993 (890) <br> response to lipid <br> 4.250e-02","GO:2000351 (53) <br> regulation of endothelial cell apoptotic process <br> 4.519e-02","GO:0031324 (2376) <br> negative regulation of cellular metabolic process <br> 4.651e-02","GO:0006109 (180) <br> regulation of carbohydrate metabolic process <br> 4.654e-02","GO:0070887 (2821) <br> cellular response to chemical stimulus <br> 4.713e-02","GO:0040011 (1646) <br> locomotion <br> 4.737e-02","GO:0051240 (1747) <br> positive regulation of multicellular organismal process <br> 4.749e-02"],"key":["GO:0065007","GO:0032501","GO:0008150","GO:0048519","GO:0065008","GO:0048523","GO:0050789","GO:0042730","GO:0071704","GO:0080090","GO:0009987","GO:0050794","GO:0060255","GO:0030195","GO:0050819","GO:1900047","GO:0019222","GO:0008152","GO:0030193","GO:1900046","GO:0050818","GO:0031323","GO:0048518","GO:0051171","GO:0042221","GO:0061045","GO:0032502","GO:0010605","GO:0044238","GO:0006807","GO:0044237","GO:0003008","GO:1902042","GO:0009892","GO:0050896","GO:1903035","GO:0048856","GO:0051179","GO:0043170","GO:2000352","GO:0051235","GO:0034116","GO:0051246","GO:1901564","GO:0007275","GO:0010883","GO:0019915","GO:1901575","GO:0010604","GO:0044262","GO:0061041","GO:0042592","GO:0009893","GO:0032879","GO:0051172","GO:0010675","GO:0031639","GO:0010033","GO:1901700","GO:0048583","GO:0051239","GO:0009894","GO:0010467","GO:1902041","GO:1904036","GO:0007596","GO:0010876","GO:0007599","GO:0050817","GO:1905952","GO:0010628","GO:0010468","GO:1903034","GO:0009895","GO:0019538","GO:0010884","GO:0048584","GO:0043086","GO:0072378","GO:1902373","GO:0034114","GO:1900048","GO:0030194","GO:0048513","GO:0033993","GO:2000351","GO:0031324","GO:0006109","GO:0070887","GO:0040011","GO:0051240"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[22.1811149141398,21.6302399060923,22.672659520516,21.3655328342547,20.99870194206,21.2646656923826,22.1206667988118,13.2269677575488,22.1127118954382,21.4420728051743,22.5189392828654,22.0729167937244,21.5067951189774,14.849607111858,14.8897757234991,14.8897757234991,21.5907472755744,22.1657486477853,15.8522097667904,15.8728459953262,15.8931913641335,21.4818775008735,21.5118291424317,21.4118093221725,21.0879119004744,15.5532960729524,21.5112263303276,20.6159088316184,22.0486579001804,21.9882346936778,22.0658759487636,20.0177805743116,14.0082800347689,20.7148249558513,21.8258163472951,15.9132536012379,21.4243246011495,21.507198728705,21.9486141444188,14.1411225406242,17.9643439044818,12.5628666298872,20.645541325712,21.5677649508492,21.328194169814,15.3337012846078,16.0822112162691,20.2268440099866,20.8910800984456,17.8212275864661,16.7246403836453,20.1724338444707,20.9898856212222,20.6459910287057,20.4448119773284,16.8525844149833,13.3339865625059,20.7675813675783,20.0866651799305,20.9879166213135,20.7807794493014,19.2589235236653,21.3897137252826,14.928921455919,14.928921455919,16.9880534913844,18.2590735700572,17.0063430455973,17.0063430455973,17.0153954918685,20.3906477904016,21.1523888526142,17.077110792071,17.8015380021475,21.4084910419731,13.7862441918272,20.345571017088,19.0501990078366,11.5015907162123,15.1448099527621,13.8636394144686,13.8636394144686,13.8636394144686,20.8637994591789,19.3282464562775,15.1779846322092,20.533419057767,17.136109390795,20.7353225700316,20.0933286855671,20.1655303481851],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[64.0726017333845,50.8928555223487,51.0306734519082,51.0524341776282,51.4259933024869,50.8964823099687,57.7583644869857,57.7692448498456,50.9001090975887,57.7619912746056,61.519343248913,59.7821119789389,57.7764984250856],"y":[12.6654552456568,11.1503059453345,4.62476001696794,4.19651831969768,4.17015407570673,3.76629110764122,3.33431854758353,2.93681533396207,2.70050202491132,2.64696342507917,2.4595767889025,2.35218317443729,1.49960633200339],"text":["GO:0110165 (15331) <br> cellular anatomical entity <br> 2.160e-13","GO:0005575 (16333) <br> cellular_component <br> 7.074e-12","GO:0005615 (1358) <br> extracellular space <br> 2.373e-05","GO:0005622 (12409) <br> intracellular <br> 6.360e-05","GO:0005737 (9997) <br> cytoplasm <br> 6.758e-05","GO:0005576 (2069) <br> extracellular region <br> 1.713e-04","GO:0043226 (11287) <br> organelle <br> 4.631e-04","GO:0043229 (10985) <br> intracellular organelle <br> 1.157e-03","GO:0005577 (6) <br> fibrinogen complex <br> 1.993e-03","GO:0043227 (10183) <br> membrane-bounded organelle <br> 2.254e-03","GO:0072562 (7) <br> blood microparticle <br> 3.471e-03","GO:0062023 (329) <br> collagen-containing extracellular matrix <br> 4.444e-03","GO:0043231 (9530) <br> intracellular membrane-bounded organelle <br> 3.165e-02"],"key":["GO:0110165","GO:0005575","GO:0005615","GO:0005622","GO:0005737","GO:0005576","GO:0043226","GO:0043229","GO:0005577","GO:0043227","GO:0072562","GO:0062023","GO:0043231"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(16,150,24,1)","opacity":0.8,"size":[22.6106378719748,22.6771653543307,19.857924281073,22.3867151394929,22.1550144889914,20.3689149305258,22.2854776116061,22.2564038814659,10.3791711622863,22.174890394537,10.8269890351464,18.0030523997013,22.1033332412498],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(16,150,24,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"GO:CC","legendgroup":"GO:CC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[8.81376447513282,8.87541109358799,3.4686400279025,8.82464329015432,9.60429170002849],"y":[8.03945730259027,7.5235508310413,6.894264659729,3.73327970277724,2.50469775922711],"text":["GO:0005488 (12742) <br> binding <br> 9.132e-09","GO:0005515 (9230) <br> protein binding <br> 2.995e-08","GO:0003674 (16250) <br> molecular_function <br> 1.276e-07","GO:0005496 (101) <br> steroid binding <br> 1.848e-04","GO:0008289 (723) <br> lipid binding <br> 3.128e-03"],"key":["GO:0005488","GO:0005515","GO:0003674","GO:0005496","GO:0008289"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[22.4149055418028,22.0686981033598,22.6718204691843,16.2490716540754,19.061012747482],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[185.182638165491,184.107350840567],"y":[7.74042647632453,2.05342945977135],"text":["KEGG:04610 (85) <br> Complement and coagulation cascades <br> 1.818e-08","KEGG:00000 (6724) <br> KEGG root term <br> 8.842e-03"],"key":["KEGG:04610","KEGG:00000"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(221,68,119,1)","opacity":0.8,"size":[15.9718144123979,21.722073290074],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(221,68,119,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"KEGG","legendgroup":"KEGG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[190.758040565264,193.983488491609,191.447393890355,193.446518533118,193.943578562262,195.329541563234,192.064183707541,195.155389144263,193.377583200609,192.343553212972,193.40660860377,193.395724077585,194.154012735184,195.60528289327,191.657828063277,193.39935225298,191.900915814756,194.161269085975,193.373955025213],"y":[4.32897922957489,3.31494810647701,3.24315901652088,2.98391423816892,2.81005491126212,2.53606766039203,2.42308453649021,2.30331672820784,2.02049136741128,1.91867812074311,1.70621156120837,1.69373634856089,1.61582108983224,1.56646031955938,1.56646031955938,1.37846519032636,1.36788394809692,1.32889562418354,1.30137451009361],"text":["REAC:R-MMU-140875 (20) <br> Common Pathway of Fibrin Clot Formation <br> 4.688e-05","REAC:R-MMU-5686938 (14) <br> Regulation of TLR by endogenous ligand <br> 4.842e-04","REAC:R-MMU-140877 (32) <br> Formation of Fibrin Clot (Clotting Cascade) <br> 5.713e-04","REAC:R-MMU-8957275 (97) <br> Post-translational protein phosphorylation <br> 1.038e-03","REAC:R-MMU-381426 (103) <br> Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs) <br> 1.549e-03","REAC:R-MMU-382551 (598) <br> Transport of small molecules <br> 2.910e-03","REAC:R-MMU-168249 (851) <br> Innate Immune System <br> 3.775e-03","REAC:R-MMU-168898 (123) <br> Toll-like Receptor Cascades <br> 4.974e-03","REAC:R-MMU-174824 (56) <br> Plasma lipoprotein assembly, remodeling, and clearance <br> 9.539e-03","REAC:R-MMU-8964041 (2) <br> LDL remodeling <br> 1.206e-02","REAC:R-MMU-114608 (105) <br> Platelet degranulation  <br> 1.967e-02","REAC:R-MMU-76009 (34) <br> Platelet Aggregation (Plug Formation) <br> 2.024e-02","REAC:R-MMU-76005 (109) <br> Response to elevated platelet cytosolic Ca2+ <br> 2.422e-02","REAC:R-MMU-372708 (14) <br> p130Cas linkage to MAPK signaling for integrins <br> 2.714e-02","REAC:R-MMU-354194 (14) <br> GRB2:SOS provides linkage to MAPK signaling for Integrins  <br> 2.714e-02","REAC:R-MMU-76002 (230) <br> Platelet activation, signaling and aggregation <br> 4.183e-02","REAC:R-MMU-109582 (511) <br> Hemostasis <br> 4.287e-02","REAC:R-MMU-975634 (42) <br> Retinoid metabolism and transport <br> 4.689e-02","REAC:R-MMU-8963898 (17) <br> Plasma lipoprotein assembly <br> 4.996e-02"],"key":["REAC:R-MMU-140875","REAC:R-MMU-5686938","REAC:R-MMU-140877","REAC:R-MMU-8957275","REAC:R-MMU-381426","REAC:R-MMU-382551","REAC:R-MMU-168249","REAC:R-MMU-168898","REAC:R-MMU-174824","REAC:R-MMU-8964041","REAC:R-MMU-114608","REAC:R-MMU-76009","REAC:R-MMU-76005","REAC:R-MMU-372708","REAC:R-MMU-354194","REAC:R-MMU-76002","REAC:R-MMU-109582","REAC:R-MMU-975634","REAC:R-MMU-8963898"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(51,102,204,1)","opacity":0.8,"size":[13.3339865625059,12.5628666298872,14.2638648249419,16.184666169372,16.2802044315257,18.8127672691082,19.2710142596048,16.558505271697,15.2733393369483,3.77952755905512,16.3106633343576,14.3778687414291,16.3696673575193,12.5628666298872,12.5628666298872,17.4950785556368,18.6040055265771,14.7659892698552,12.9906351928736],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(51,102,204,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"REAC","legendgroup":"REAC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[199.726645872062,199.289308867157],"y":[6.97610538918962,5.80813620292009],"text":["WP:WP460 (20) <br> Blood Clotting Cascade <br> 1.057e-07","WP:000000 (3972) <br> WIKIPATHWAYS <br> 1.555e-06"],"key":["WP:WP460","WP:000000"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,153,198,1)","opacity":0.8,"size":[13.3339865625059,21.1307493164437],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,153,198,1)"}},"hoveron":"points","set":"SharedData5e8be6f1","name":"WP","legendgroup":"WP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[2,46.3456762993078],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[49.9716514668718,66.1580046148775],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[69.7839797824415,180.481375673003],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[184.107350840567,186.058125480716],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(221,68,119,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[189.68410064828,195.663333699593],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[199.289308867157,200.003625975168],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(0,153,198,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":71.8183902982224,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"tickvals":[24.1728381496539,58.0648280408746,125.132677727722,185.082738160642,192.673717173937,199.646467421163],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0217413864674139,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"7a2318fd2491":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"7a23615fd696":{"x":{},"xend":{},"y":{},"yend":{}},"7a237a27ca8":{"x":{},"xend":{},"y":{},"yend":{}},"7a23966ed9e":{"x":{},"xend":{},"y":{},"yend":{}},"7a2315551896":{"x":{},"xend":{},"y":{},"yend":{}},"7a2321e49c3d":{"x":{},"xend":{},"y":{},"yend":{}},"7a237d50bd50":{"x":{},"xend":{},"y":{},"yend":{}},"7a235caf850d":{"x":{},"y":{}},"7a2367a8d945":{"yintercept":{}}},"cur_data":"7a2318fd2491","visdat":{"7a2318fd2491":["function (y) ","x"],"7a23615fd696":["function (y) ","x"],"7a237a27ca8":["function (y) ","x"],"7a23966ed9e":["function (y) ","x"],"7a2315551896":["function (y) ","x"],"7a2321e49c3d":["function (y) ","x"],"7a237d50bd50":["function (y) ","x"],"7a235caf850d":["function (y) ","x"],"7a2367a8d945":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedData5e8be6f1"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

### pattern3.FDR0.05 (Lingon specific)


```r
gost_pattern3 <- ggo(pattern3.FDR0.05)
gostplot(gost_pattern3)
```

<!--html_preserve--><div id="htmlwidget-ec3abe2be1024ad27d2f" style="width:960px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-ec3abe2be1024ad27d2f">{"x":{"data":[{"x":[112.724184254927,79.1683109068117,116.796287752866,114.917971090415],"y":[1.9584746531904,1.8407288235717,1.8407288235717,1.64914863120182],"text":["GO:0042438 (22) <br> melanin biosynthetic process <br> 1.100e-02","GO:0006582 (24) <br> melanin metabolic process <br> 1.443e-02","GO:0044550 (24) <br> secondary metabolite biosynthetic process <br> 1.443e-02","GO:0043473 (88) <br> pigmentation <br> 2.243e-02"],"key":["GO:0042438","GO:0006582","GO:0044550","GO:0043473"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[13.5297239329324,13.705048182282,13.705048182282,16.0280782195047],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedDataf7a43ec5","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[50.8928555223487,64.0726017333845,55.9014492255519,58.835520410122,61.8167398337521],"y":[2.23818328053478,1.98542684680461,1.308035390324,1.308035390324,1.308035390324],"text":["GO:0005575 (16333) <br> cellular_component <br> 5.779e-03","GO:0110165 (15331) <br> cellular anatomical entity <br> 1.034e-02","GO:0033162 (14) <br> melanosome membrane <br> 4.920e-02","GO:0045009 (14) <br> chitosome <br> 4.920e-02","GO:0090741 (14) <br> pigment granule membrane <br> 4.920e-02"],"key":["GO:0005575","GO:0110165","GO:0033162","GO:0045009","GO:0090741"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(16,150,24,1)","opacity":0.8,"size":[22.6771653543307,22.6106378719748,12.5628666298872,12.5628666298872,12.5628666298872],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(16,150,24,1)"}},"hoveron":"points","set":"SharedDataf7a43ec5","name":"GO:CC","legendgroup":"GO:CC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[6.03966664465033,15.4861043549861],"y":[3.04068431063149,2.56399597296634],"text":["GO:0004503 (2) <br> monophenol monooxygenase activity <br> 9.106e-04","GO:0016716 (3) <br> oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, another compound as one donor, and incorporation of one atom of oxygen <br> 2.729e-03"],"key":["GO:0004503","GO:0016716"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[3.77952755905512,7.78888707317388],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedDataf7a43ec5","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[185.367907265394],"y":[1.57577915186593],"text":"KEGG:04916 (95) <br> Melanogenesis <br> 2.656e-02","key":["KEGG:04916"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(221,68,119,1)","opacity":0.8,"size":16.1513294155561,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(221,68,119,1)"}},"hoveron":"points","set":"SharedDataf7a43ec5","name":"KEGG","legendgroup":"KEGG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isSimpleKey":true,"_isNestedKey":false,"frame":null},{"x":[192.506821105757],"y":[2.2573008263587],"text":"REAC:R-MMU-5662702 (5) <br> Melanin biosynthesis <br> 5.530e-03","key":["REAC:R-MMU-5662702"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(51,102,204,1)","opacity":0.8,"size":9.80671785770929,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(51,102,204,1)"}},"hoveron":"points","set":"SharedDataf7a43ec5","name":"REAC","legendgroup":"REAC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isSimpleKey":true,"_isNestedKey":false,"frame":null},{"x":[2,46.3456762993078],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[49.9716514668718,66.1580046148775],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[69.7839797824415,180.481375673003],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[184.107350840567,186.058125480716],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(221,68,119,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[189.68410064828,195.663333699593],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[199.289308867157,200.003625975168],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(0,153,198,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":71.8183902982224,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"tickvals":[24.1728381496539,58.0648280408746,125.132677727722,185.082738160642,192.673717173937,199.646467421163],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0217413864674139,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"7a23e8bf82b":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"7a2315515ab4":{"x":{},"xend":{},"y":{},"yend":{}},"7a236eedff4":{"x":{},"xend":{},"y":{},"yend":{}},"7a231e814ef6":{"x":{},"xend":{},"y":{},"yend":{}},"7a2331d31104":{"x":{},"xend":{},"y":{},"yend":{}},"7a2325c7951f":{"x":{},"xend":{},"y":{},"yend":{}},"7a2354a5a771":{"x":{},"xend":{},"y":{},"yend":{}},"7a237d353638":{"x":{},"y":{}},"7a236f38090e":{"yintercept":{}}},"cur_data":"7a23e8bf82b","visdat":{"7a23e8bf82b":["function (y) ","x"],"7a2315515ab4":["function (y) ","x"],"7a236eedff4":["function (y) ","x"],"7a231e814ef6":["function (y) ","x"],"7a2331d31104":["function (y) ","x"],"7a2325c7951f":["function (y) ","x"],"7a2354a5a771":["function (y) ","x"],"7a237d353638":["function (y) ","x"],"7a236f38090e":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedDataf7a43ec5"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

### pattern4.FDR0.05 (HFD specific)


```r
gost_pattern4 <- ggo(pattern4.FDR0.05)
gostplot(gost_pattern4)
```

<!--html_preserve--><div id="htmlwidget-e82c15034d24ce25fe01" style="width:960px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-e82c15034d24ce25fe01">{"x":{"data":[{"x":[128.414292746097,125.219703982353,111.618225602325,80.296026123071,101.947433056461,82.7073785951368,85.4233229584108,76.8077237499473,91.0546468518215,101.951059150404,125.172564761095,107.350313031466,143.7019048096,127.286577529838,80.5571048869638,106.752007530878,124.748311769769,150.319526255495,97.5018418823973,69.7912319703274,150.55884845573,132.225317480144,132.196308728601,71.4229742446576,117.115384019846,81.3874803999008,127.029124859888,128.501319000728,85.2057573218335,87.1529697692009,99.9567074817786,71.5426353447752,71.4157220567717,125.847018234485,124.03759735695,125.832513858713,87.2472482117177,80.9849839722327,147.208337652439,97.5489811036558,127.00736829623,147.201085464553,124.055727826664,137.842136997783,117.267679965451,70.1683457403948,143.310286663761,142.693850693458],"y":[7.90286571056473,7.90286571056473,6.75389118241118,6.68253698690871,5.17754853090855,4.8513013643868,4.04280832469327,4.00408341918802,3.9714770800843,3.64933808736568,3.45127696344242,3.34412447976156,3.29753689161745,3.00590792221126,2.97006940151466,2.94095633527592,2.84318269861266,2.81462999237966,2.72358233103101,2.71822917840963,2.64909820838986,2.6271507234876,2.46741025667015,2.46741025667015,2.11190792040779,2.01027334287912,2.00625715225481,1.99463248404688,1.96661773647802,1.95402534257319,1.87754685561566,1.87197082997098,1.77720857462576,1.75896987573846,1.7436511684503,1.70911778249309,1.68651200610527,1.68069983963363,1.66120562587646,1.63062724155544,1.59548623880262,1.5865513145401,1.56731116123891,1.41407096965972,1.38473816331537,1.3261984887232,1.30955882017319,1.30955882017319],"text":["GO:0051674 (1483) <br> localization of cell <br> 1.251e-08","GO:0048870 (1483) <br> cell motility <br> 1.251e-08","GO:0040011 (1646) <br> locomotion <br> 1.762e-07","GO:0006928 (1872) <br> movement of cell or subcellular component <br> 2.077e-07","GO:0032501 (6189) <br> multicellular organismal process <br> 6.644e-06","GO:0008150 (16263) <br> biological_process <br> 1.408e-05","GO:0009987 (14055) <br> cellular process <br> 9.061e-05","GO:0003341 (144) <br> cilium movement <br> 9.906e-05","GO:0016477 (1348) <br> cell migration <br> 1.068e-04","GO:0032502 (5562) <br> developmental process <br> 2.242e-04","GO:0048856 (5147) <br> anatomical structure development <br> 3.538e-04","GO:0035295 (1035) <br> tube development <br> 4.528e-04","GO:0072359 (1042) <br> circulatory system development <br> 5.040e-04","GO:0051270 (986) <br> regulation of cellular component movement <br> 9.865e-04","GO:0007017 (798) <br> microtubule-based process <br> 1.071e-03","GO:0035082 (60) <br> axoneme assembly <br> 1.146e-03","GO:0048731 (4264) <br> system development <br> 1.435e-03","GO:0098609 (726) <br> cell-cell adhesion <br> 1.532e-03","GO:0022414 (1350) <br> reproductive process <br> 1.890e-03","GO:0000003 (1351) <br> reproduction <br> 1.913e-03","GO:0098742 (208) <br> cell-cell adhesion via plasma-membrane adhesion molecules <br> 2.243e-03","GO:0060294 (107) <br> cilium movement involved in cell motility <br> 2.360e-03","GO:0060285 (113) <br> cilium-dependent cell motility <br> 3.409e-03","GO:0001539 (113) <br> cilium or flagellum-dependent cell motility <br> 3.409e-03","GO:0044703 (923) <br> multi-organism reproductive process <br> 7.728e-03","GO:0007275 (4726) <br> multicellular organism development <br> 9.766e-03","GO:0051179 (5542) <br> localization <br> 9.857e-03","GO:0051704 (942) <br> multi-organism process <br> 1.012e-02","GO:0009893 (3513) <br> positive regulation of metabolic process <br> 1.080e-02","GO:0010604 (3225) <br> positive regulation of macromolecule metabolic process <br> 1.112e-02","GO:0031325 (3106) <br> positive regulation of cellular metabolic process <br> 1.326e-02","GO:0001578 (91) <br> microtubule bundle formation <br> 1.343e-02","GO:0001525 (498) <br> angiogenesis <br> 1.670e-02","GO:0050794 (9266) <br> regulation of cellular process <br> 1.742e-02","GO:0048513 (3150) <br> animal organ development <br> 1.804e-02","GO:0050789 (9684) <br> regulation of biological process <br> 1.954e-02","GO:0010631 (273) <br> epithelial cell migration <br> 2.058e-02","GO:0007155 (1215) <br> cell adhesion <br> 2.086e-02","GO:0090132 (275) <br> epithelium migration <br> 2.182e-02","GO:0022610 (1225) <br> biological adhesion <br> 2.341e-02","GO:0051173 (2912) <br> positive regulation of nitrogen compound metabolic process <br> 2.538e-02","GO:0090130 (281) <br> tissue migration <br> 2.591e-02","GO:0048518 (5565) <br> positive regulation of biological process <br> 2.708e-02","GO:0065007 (10242) <br> biological regulation <br> 3.854e-02","GO:0044782 (298) <br> cilium organization <br> 4.123e-02","GO:0000226 (553) <br> microtubule cytoskeleton organization <br> 4.718e-02","GO:0072240 (2) <br> metanephric DCT cell differentiation <br> 4.903e-02","GO:0072069 (2) <br> DCT cell differentiation <br> 4.903e-02"],"key":["GO:0051674","GO:0048870","GO:0040011","GO:0006928","GO:0032501","GO:0008150","GO:0009987","GO:0003341","GO:0016477","GO:0032502","GO:0048856","GO:0035295","GO:0072359","GO:0051270","GO:0007017","GO:0035082","GO:0048731","GO:0098609","GO:0022414","GO:0000003","GO:0098742","GO:0060294","GO:0060285","GO:0001539","GO:0044703","GO:0007275","GO:0051179","GO:0051704","GO:0009893","GO:0010604","GO:0031325","GO:0001578","GO:0001525","GO:0050794","GO:0048513","GO:0050789","GO:0010631","GO:0007155","GO:0090132","GO:0022610","GO:0051173","GO:0090130","GO:0048518","GO:0065007","GO:0044782","GO:0000226","GO:0072240","GO:0072069"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(255,153,0,1)","opacity":0.8,"size":[19.9661210354986,19.9661210354986,20.0933286855671,20.2489203596558,21.6302399060923,22.672659520516,22.5189392828654,16.8007251732528,19.8488093985205,21.5112263303276,21.4243246011495,19.5194957465375,19.5279824452683,19.458295568636,19.1885119254307,15.3917186760792,21.2116038419213,19.066383166014,19.8506381869434,19.8515514874927,17.3489898859544,16.3404757245714,16.4262840621335,16.4262840621335,19.3745932001494,21.328194169814,21.507198728705,19.4004716889662,20.9898856212222,20.8910800984456,20.8474700399405,16.0822112162691,18.5695071271461,22.0729167937244,20.8637994591789,22.1206667988118,17.7405920191581,19.7201518048161,17.7509522214783,19.7303415895367,20.7723991542435,17.7815418334201,21.5118291424317,22.1811149141398,17.8644538114196,18.7092529725439,3.77952755905512,3.77952755905512],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(255,153,0,1)"}},"hoveron":"points","set":"SharedData0f94fd85","name":"GO:BP","legendgroup":"GO:BP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[64.0871088838644,57.6350537079061,50.8928555223487,64.0726017333845,52.0389204102648,61.827620196612,51.8140595778256,52.0352936226449,63.8296069628453,55.7019759064526,53.7580177421391,57.7801252127056,57.7656180622256,58.9044293749018,53.1632245724611,53.9901321498183,53.0725548819614,62.7850921282888,63.3762585103468,61.4431807088933,54.9947523205549,51.8974756930853,63.0824887131278],"y":[9.061337501164,8.52935537201078,7.16030606862717,7.1331220624541,5.83197717899357,5.76183724122391,5.03278915268572,4.00889979367871,3.74046359137686,3.45083206018332,2.86370588699054,2.79459751184077,2.75579024028252,2.4290481617979,2.42898777021701,2.42807515034065,2.03512567156358,2.01380875742482,2.00278874547579,2.00273954845783,1.83305280787368,1.54399002310999,1.51409349867071],"text":["GO:0120025 (2163) <br> plasma membrane bounded cell projection <br> 8.683e-10","GO:0042995 (2252) <br> cell projection <br> 2.956e-09","GO:0005575 (16333) <br> cellular_component <br> 6.913e-08","GO:0110165 (15331) <br> cellular anatomical entity <br> 7.360e-08","GO:0005930 (109) <br> axoneme <br> 1.472e-06","GO:0097014 (111) <br> ciliary plasm <br> 1.730e-06","GO:0005856 (1968) <br> cytoskeleton <br> 9.273e-06","GO:0005929 (589) <br> cilium <br> 9.797e-05","GO:0099568 (248) <br> cytoplasmic region <br> 1.818e-04","GO:0032838 (206) <br> plasma membrane bounded cell projection cytoplasm <br> 3.541e-04","GO:0030054 (1836) <br> cell junction <br> 1.369e-03","GO:0043232 (3941) <br> intracellular non-membrane-bounded organelle <br> 1.605e-03","GO:0043228 (3955) <br> non-membrane-bounded organelle <br> 1.755e-03","GO:0045177 (427) <br> apical part of cell <br> 3.724e-03","GO:0016324 (348) <br> apical plasma membrane <br> 3.724e-03","GO:0030286 (58) <br> dynein complex <br> 3.732e-03","GO:0015630 (1185) <br> microtubule cytoskeleton <br> 9.223e-03","GO:0097729 (117) <br> 9+2 motile cilium <br> 9.687e-03","GO:0098862 (174) <br> cluster of actin-based cell projections <br> 9.936e-03","GO:0071944 (4403) <br> cell periphery <br> 9.937e-03","GO:0031514 (185) <br> motile cilium <br> 1.469e-02","GO:0005886 (4281) <br> plasma membrane <br> 2.858e-02","GO:0098590 (1181) <br> plasma membrane region <br> 3.061e-02"],"key":["GO:0120025","GO:0042995","GO:0005575","GO:0110165","GO:0005930","GO:0097014","GO:0005856","GO:0005929","GO:0099568","GO:0032838","GO:0030054","GO:0043232","GO:0043228","GO:0045177","GO:0016324","GO:0030286","GO:0015630","GO:0097729","GO:0098862","GO:0071944","GO:0031514","GO:0005886","GO:0098590"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(16,150,24,1)","opacity":0.8,"size":[20.4219211126518,20.4698804912106,22.6771653543307,22.6106378719748,16.3696673575193,16.3982625041893,20.3090042381597,18.7927579106113,17.6035509689458,17.3348679917457,20.225531780848,21.1217956425754,21.1258484990633,18.3618917074693,18.0810854513317,15.3337012846078,19.6890312047081,16.4806917307213,17.0857000429034,21.2480437232029,17.1767112481104,21.2161279031282,19.684817704093],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(16,150,24,1)"}},"hoveron":"points","set":"SharedData0f94fd85","name":"GO:CC","legendgroup":"GO:CC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[3.4686400279025,10.2969095897307,8.81376447513282,8.87541109358799,35.7859731851051,28.3629950687683,46.2405144207666],"y":[6.17606734763606,2.92837326807398,2.34090675271073,2.30001298375351,2.11985878616159,2.05709838346046,1.82763278307357],"text":["GO:0003674 (16250) <br> molecular_function <br> 6.667e-07","GO:0008569 (18) <br> ATP-dependent microtubule motor activity, minus-end-directed <br> 1.179e-03","GO:0005488 (12742) <br> binding <br> 4.561e-03","GO:0005515 (9230) <br> protein binding <br> 5.012e-03","GO:0051959 (28) <br> dynein light intermediate chain binding <br> 7.588e-03","GO:0045505 (29) <br> dynein intermediate chain binding <br> 8.768e-03","GO:1990939 (33) <br> ATP-dependent microtubule motor activity <br> 1.487e-02"],"key":["GO:0003674","GO:0008569","GO:0005488","GO:0005515","GO:0051959","GO:0045505","GO:1990939"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(220,57,18,1)","opacity":0.8,"size":[22.6718204691843,13.1128330495735,22.4149055418028,22.0686981033598,14.0082800347689,14.0760608836278,14.3218846504567],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(220,57,18,1)"}},"hoveron":"points","set":"SharedData0f94fd85","name":"GO:MF","legendgroup":"GO:MF","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[184.924714516607,184.939245426404,185.364274537944,185.360641810495,185.139045436102,185.411499994782,185.647627278972,184.997369065589,185.698485463259,185.683954553462,184.107350840567,185.164474528246],"y":[4.28044363521996,3.32312642475726,2.70556108099073,2.37013520487837,2.35077918394624,2.01871274594028,1.80700790032793,1.77869203668572,1.66010122525358,1.43797504043882,1.40510803438361,1.30394534558257],"text":["KEGG:04015 (196) <br> Rap1 signaling pathway <br> 5.243e-05","KEGG:04024 (193) <br> cAMP signaling pathway <br> 4.752e-04","KEGG:04915 (115) <br> Estrogen signaling pathway <br> 1.970e-03","KEGG:04914 (81) <br> Progesterone-mediated oocyte maturation <br> 4.264e-03","KEGG:04371 (133) <br> Apelin signaling pathway <br> 4.459e-03","KEGG:04928 (96) <br> Parathyroid hormone synthesis, secretion and action <br> 9.578e-03","KEGG:05200 (494) <br> Pathways in cancer <br> 1.560e-02","KEGG:04114 (108) <br> Oocyte meiosis <br> 1.665e-02","KEGG:05218 (63) <br> Melanoma <br> 2.187e-02","KEGG:05214 (72) <br> Glioma <br> 3.648e-02","KEGG:00000 (6724) <br> KEGG root term <br> 3.935e-02","KEGG:04514 (137) <br> Cell adhesion molecules <br> 4.967e-02"],"key":["KEGG:04015","KEGG:04024","KEGG:04915","KEGG:04914","KEGG:04371","KEGG:04928","KEGG:05200","KEGG:04114","KEGG:05218","KEGG:05214","KEGG:00000","KEGG:04514"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(221,68,119,1)","opacity":0.8,"size":[17.2619020696014,17.2392046697092,16.4537536562816,15.8931913641335,16.6791853660405,16.1680958215515,18.5586943215225,16.3551475853292,15.4747105262571,15.6988922791974,21.722073290074,16.7246403836453],"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(221,68,119,1)"}},"hoveron":"points","set":"SharedData0f94fd85","name":"KEGG","legendgroup":"KEGG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isNestedKey":false,"frame":null},{"x":[199.828691173206],"y":[1.94953637868267],"text":"WP:WP3654 (26) <br> Novel Jun-Dmp1 Pathway <br> 1.123e-02","key":["WP:WP3654"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(0,153,198,1)","opacity":0.8,"size":13.8636394144686,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(0,153,198,1)"}},"hoveron":"points","set":"SharedData0f94fd85","name":"WP","legendgroup":"WP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","_isSimpleKey":true,"_isNestedKey":false,"frame":null},{"x":[2,46.3456762993078],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(220,57,18,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[49.9716514668718,66.1580046148775],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(16,150,24,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[69.7839797824415,180.481375673003],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(255,153,0,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[184.107350840567,186.058125480716],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(221,68,119,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[189.68410064828,195.663333699593],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(51,102,204,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[199.289308867157,200.003625975168],"y":[-1,-1],"text":"","type":"scatter","mode":"lines","line":{"width":11.3385826771654,"color":"rgba(0,153,198,1)","dash":"solid"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[180],"y":[16.2],"text":"values above this threshold are capped","hovertext":"","textfont":{"size":7.55905511811024,"color":"rgba(190,190,190,1)"},"type":"scatter","mode":"text","hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0,210],"y":[16,16],"text":"","type":"scatter","mode":"lines","line":{"width":0.755905511811024,"color":"rgba(190,190,190,1)","dash":"dash"},"hoveron":"points","showlegend":false,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":29.2835201328352,"r":6.6417600664176,"b":71.8183902982224,"l":61.0377750103778},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0,210],"tickmode":"array","ticktext":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"tickvals":[24.1728381496539,58.0648280408746,125.132677727722,185.082738160642,192.673717173937,199.646467421163],"categoryorder":"array","categoryarray":["GO:MF","GO:CC","GO:BP","KEGG","REAC","WP"],"nticks":null,"ticks":"","tickcolor":null,"ticklen":3.65296803652968,"tickwidth":0,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":10.6268161062682},"tickangle":-45,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.132835201328352,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":"","hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-1,18.5],"tickmode":"array","ticktext":["0","2","4","6","8","10","12","14",">16"],"tickvals":[0,2,4,6,8,10,12,14,16],"categoryorder":"array","categoryarray":["0","2","4","6","8","10","12","14",">16"],"nticks":null,"ticks":"outside","tickcolor":"rgba(190,190,190,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":true,"linecolor":"rgba(190,190,190,1)","linewidth":0.66417600664176,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":"","hoverformat":".2f"},"annotations":[{"text":"-log10(p-adj)","x":-0.0217413864674139,"y":0.5,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-90,"xanchor":"right","yanchor":"center","annotationType":"axis"},{"text":"query_1","x":0.5,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(169,169,169,1)","family":"","size":13.2835201328352},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"center","yanchor":"bottom"}],"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1},{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":24.9730178497302,"yanchor":1,"ysizemode":"pixel"}],"showlegend":false,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895}},"hovermode":"closest","barmode":"relative","dragmode":"zoom"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"7a234e15135c":{"colour":{},"size":{},"alpha":{},"key":{},"x":{},"y":{},"text":{},"type":"scatter"},"7a23165f6dd0":{"x":{},"xend":{},"y":{},"yend":{}},"7a231683f37a":{"x":{},"xend":{},"y":{},"yend":{}},"7a234b7d8be5":{"x":{},"xend":{},"y":{},"yend":{}},"7a2351a4d65b":{"x":{},"xend":{},"y":{},"yend":{}},"7a234c247cf8":{"x":{},"xend":{},"y":{},"yend":{}},"7a236cdc736d":{"x":{},"xend":{},"y":{},"yend":{}},"7a23620ce31f":{"x":{},"y":{}},"7a233e593f7c":{"yintercept":{}}},"cur_data":"7a234e15135c","visdat":{"7a234e15135c":["function (y) ","x"],"7a23165f6dd0":["function (y) ","x"],"7a231683f37a":["function (y) ","x"],"7a234b7d8be5":["function (y) ","x"],"7a2351a4d65b":["function (y) ","x"],"7a234c247cf8":["function (y) ","x"],"7a236cdc736d":["function (y) ","x"],"7a23620ce31f":["function (y) ","x"],"7a233e593f7c":["function (y) ","x"]},"highlight":{"on":"plotly_click","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0,"ctGroups":["SharedData0f94fd85"]},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## Results tables

- `results/clusterprofiler/` - clusterProfiler GO and KEGG enrichment results.
- `results/gprofiler/` - gProfiler GO, KEGG, Reactome, and WikiPathways enrichment results.
