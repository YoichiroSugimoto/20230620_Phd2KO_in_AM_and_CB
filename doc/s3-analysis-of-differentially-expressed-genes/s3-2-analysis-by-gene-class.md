---
title: "s3-2 Analysis of genes associated with acute oxygen sensing"
author: "Yoichiro Sugimoto"
date: "23 June, 2023"
vignette: >
  %\VignetteIndexEntry{Bioconductor style for PDF documents}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
output:
   html_document:
     highlight: haddock
     toc: yes
     toc_depth: 2
     keep_md: yes
     fig_width: 5
     fig_height: 5
---

# Overview

Changes in the expression of genes that are associated with acute oxygen sensing capability by previous studies will be analysed.


# Setup



```r
project.dir <- "/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB"

renv::restore(file.path(project.dir, "R"))
```

```
## * The library is already synchronized with the lockfile.
```

```r
library("ggrepel")

library("org.Mm.eg.db")
library("GO.db")

processors <- 7

temp <- sapply(list.files(
    file.path(project.dir, "R/functions"),
    full.names = TRUE
), source)

set.seed(1)
```

## The input and output files and the directories


```r
annot.dir <- file.path(project.dir, "annotation/")
annot.ps.dir <- file.path(annot.dir, "mm39_annotation/processed_data/")

results.dir <- file.path(project.dir, "results")
s2.dir <- file.path(results.dir, "s2")
s3.dir <- file.path(results.dir, "s3")

create.dirs(c())
```

# Potassium channel related genes


```r
res.summary.dt <- fread(file.path(s2.dir, "all-de-results.csv"))

returnCappedVal <- function(val, log2fc.cap = 7.5){
    case_when(
        abs(val) > log2fc.cap ~ sign(val) * log2fc.cap,
        TRUE ~ val
    )
}

res.summary.dt[, `:=`(
    concordant_regulation = case_when(
        tissue_specific_genes == "CB specific" & PHD2_regulated_in_AM == "PHD2KO_induced" ~
            "CB specific and Phd2KO induced in AM",
        tissue_specific_genes == "AM specific" & PHD2_regulated_in_AM == "PHD2KO_repressed" ~
            "AM specific and Phd2KO repressed in AM"
    ),
    capped_log2fc__CB_PHD2KO = returnCappedVal(log2fc__CB_PHD2KO),
    capped_log2fc__AM_PHD2KO = returnCappedVal(log2fc__AM_PHD2KO),
    capped_log2fc__Tissue = returnCappedVal(log2fc__Tissue)
)]

go.dt <- AnnotationDbi::select(
                            org.Mm.eg.db,
                            keys(org.Mm.eg.db, "GO"), "ENSEMBL", "GO"
                        ) %>%
    data.table()
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
setnames(go.dt, old = "ENSEMBL", new = "gene_id")

go.dt[, `:=`(
    gene_group = case_when(
        GO %in% c(
                    "GO:0005267"
                ) ~ "potassium channel"
    )
)]

go.dt <- go.dt[!duplicated(paste(gene_id, gene_group))]

go.dt <- go.dt[!is.na(gene_id) & !is.na(gene_group), .(gene_id, gene_group)]

collapsed.class.dt <- go.dt[
  , lapply(.SD, paste, collapse = ", "),
    by = list(gene_id), .SDcols = "gene_group"
]

res.by.class.dt <- merge(
    collapsed.class.dt,
    res.summary.dt,
    by = "gene_id",
    all.y = TRUE
)

res.by.class.dt <- res.by.class.dt[
  , tissue_vs_AM_PHD2KO_comparison_flag :=
        !is.na(log2fc__Tissue) & !is.na(padj__AM_PHD2KO)
]


## Export results

export.dt <- copy(res.by.class.dt)[
   ,.(
        gene_id, gene_name, gene_type, gene_group,
        TPM_WT__CB_PHD2KO, TPM_PHD2KO__CB_PHD2KO, TPM_WT__AM_PHD2KO, TPM_PHD2KO__AM_PHD2KO,
        log2fc__Tissue, padj__Tissue, tissue_specific_genes,
        log2fc__AM_PHD2KO, padj__AM_PHD2KO, PHD2_regulated_in_AM,
        log2fc__CB_PHD2KO, padj__CB_PHD2KO, PHD2_regulated_in_CB
    )]

setnames(
    export.dt,
    old = c(
        "TPM_WT__CB_PHD2KO", "TPM_PHD2KO__CB_PHD2KO",
        "TPM_WT__AM_PHD2KO", "TPM_PHD2KO__AM_PHD2KO",
        "gene_group"
    ),
    new = c("TPM_CB_WT", "TPM_CB_PHD2KO", "TPM_AM_WT", "TPM_AM_PHD2KO", "potassium_channel")
)

export.dt[, potassium_channel := case_when(
                potassium_channel == "potassium channel" ~ "Yes",
                TRUE ~ ""
            )]

fwrite(
    export.dt,
    file.path(s3.dir, "summary-of-differentially-expressed-genes.csv")
)
```


# Classes of genes associated with acute oxygen sensing by previous studies



```r
acute.oxygen.genes.dt <- fread(
    file.path(project.dir, "data/ref/candidate_genes_in_acute_oxygen_sensing.csv")
)



sl.dt <- export.dt[
    acute.oxygen.genes.dt[
        gene_id != "",
        list(
            gene_id = str_split_fixed(gene_id, "\\.", n = 2)[, 1],
            Pathway                                
        )]
]

fwrite(sl.dt, file.path(s3.dir, "acute_oxygen_genes.csv"))
```



# Session information


```r
sessioninfo::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R version 4.2.1 (2022-06-23)
##  os       CentOS Linux 7 (Core)
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_GB.UTF-8
##  ctype    en_GB.UTF-8
##  tz       Europe/Berlin
##  date     2023-06-23
##  pandoc   2.19.2 @ /fast/AG_Sugimoto/home/users/yoichiro/software/miniconda3/envs/20220601_CB_AM_PHD2/bin/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package          * version  date (UTC) lib source
##  AnnotationDbi    * 1.60.0   2022-11-01 [1] Bioconductor
##  Biobase          * 2.56.0   2022-04-26 [1] Bioconductor
##  BiocGenerics     * 0.42.0   2022-04-26 [1] Bioconductor
##  BiocManager        1.30.18  2022-05-18 [1] CRAN (R 4.2.1)
##  Biostrings         2.64.1   2022-08-18 [1] Bioconductor
##  bit                4.0.4    2020-08-04 [1] CRAN (R 4.2.1)
##  bit64              4.0.5    2020-08-30 [1] CRAN (R 4.2.1)
##  bitops             1.0-7    2021-04-24 [1] CRAN (R 4.2.1)
##  blob               1.2.3    2022-04-10 [1] CRAN (R 4.2.1)
##  bslib              0.4.0    2022-07-16 [1] CRAN (R 4.2.1)
##  cachem             1.0.6    2021-08-19 [1] CRAN (R 4.2.1)
##  cli                3.4.1    2022-09-23 [1] CRAN (R 4.2.1)
##  colorspace         2.0-3    2022-02-21 [1] CRAN (R 4.2.1)
##  crayon             1.5.2    2022-09-29 [1] CRAN (R 4.2.1)
##  data.table       * 1.14.4   2022-10-17 [1] CRAN (R 4.2.1)
##  DBI                1.1.3    2022-06-18 [1] CRAN (R 4.2.1)
##  digest             0.6.30   2022-10-18 [1] CRAN (R 4.2.1)
##  dplyr            * 1.0.10   2022-09-01 [1] CRAN (R 4.2.1)
##  evaluate           0.17     2022-10-07 [1] CRAN (R 4.2.1)
##  fansi              1.0.3    2022-03-24 [1] CRAN (R 4.2.1)
##  fastmap            1.1.0    2021-01-25 [1] CRAN (R 4.2.1)
##  generics           0.1.3    2022-07-05 [1] CRAN (R 4.2.1)
##  GenomeInfoDb       1.32.4   2022-09-06 [1] Bioconductor
##  GenomeInfoDbData   1.2.8    2022-10-21 [1] Bioconductor
##  ggplot2          * 3.3.6    2022-05-03 [1] CRAN (R 4.2.1)
##  ggrepel          * 0.9.1    2021-01-15 [1] CRAN (R 4.2.1)
##  glue               1.6.2    2022-02-24 [1] CRAN (R 4.2.1)
##  GO.db            * 3.16.0   2022-12-07 [1] Bioconductor
##  gtable             0.3.1    2022-09-01 [1] CRAN (R 4.2.1)
##  htmltools          0.5.3    2022-07-18 [1] CRAN (R 4.2.1)
##  httr               1.4.4    2022-08-17 [1] CRAN (R 4.2.1)
##  IRanges          * 2.30.1   2022-08-18 [1] Bioconductor
##  jquerylib          0.1.4    2021-04-26 [1] CRAN (R 4.2.1)
##  jsonlite           1.8.2    2022-10-02 [1] CRAN (R 4.2.1)
##  KEGGREST           1.36.3   2022-07-12 [1] Bioconductor
##  khroma           * 1.9.0    2022-06-18 [1] CRAN (R 4.2.1)
##  knitr            * 1.40     2022-08-24 [1] CRAN (R 4.2.1)
##  lifecycle          1.0.3    2022-10-07 [1] CRAN (R 4.2.1)
##  magrittr         * 2.0.3    2022-03-30 [1] CRAN (R 4.2.1)
##  memoise            2.0.1    2021-11-26 [1] CRAN (R 4.2.1)
##  munsell            0.5.0    2018-06-12 [1] CRAN (R 4.2.1)
##  org.Mm.eg.db     * 3.16.0   2022-12-07 [1] Bioconductor
##  pillar             1.8.1    2022-08-19 [1] CRAN (R 4.2.1)
##  pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 4.2.1)
##  png                0.1-7    2013-12-03 [1] CRAN (R 4.2.1)
##  R6                 2.5.1    2021-08-19 [1] CRAN (R 4.2.1)
##  Rcpp               1.0.9    2022-07-08 [1] CRAN (R 4.2.1)
##  RCurl              1.98-1.9 2022-10-03 [1] CRAN (R 4.2.1)
##  renv               0.16.0   2022-09-29 [1] CRAN (R 4.2.1)
##  rlang              1.0.6    2022-09-24 [1] CRAN (R 4.2.1)
##  rmarkdown        * 2.17     2022-10-07 [1] CRAN (R 4.2.1)
##  RSQLite            2.2.18   2022-10-04 [1] CRAN (R 4.2.1)
##  S4Vectors        * 0.34.0   2022-04-26 [1] Bioconductor
##  sass               0.4.2    2022-07-16 [1] CRAN (R 4.2.1)
##  scales             1.2.1    2022-08-20 [1] CRAN (R 4.2.1)
##  sessioninfo        1.2.2    2021-12-06 [1] CRAN (R 4.2.1)
##  stringi            1.7.8    2022-07-11 [1] CRAN (R 4.2.1)
##  stringr          * 1.4.1    2022-08-20 [1] CRAN (R 4.2.1)
##  tibble             3.1.8    2022-07-22 [1] CRAN (R 4.2.1)
##  tidyselect         1.2.0    2022-10-10 [1] CRAN (R 4.2.1)
##  utf8               1.2.2    2021-07-24 [1] CRAN (R 4.2.1)
##  vctrs              0.4.2    2022-09-29 [1] CRAN (R 4.2.1)
##  withr              2.5.0    2022-03-03 [1] CRAN (R 4.2.1)
##  xfun               0.34     2022-10-18 [1] CRAN (R 4.2.1)
##  XVector            0.36.0   2022-04-26 [1] Bioconductor
##  yaml               2.3.6    2022-10-18 [1] CRAN (R 4.2.1)
##  zlibbioc           1.42.0   2022-04-26 [1] Bioconductor
## 
##  [1] /fast/AG_Sugimoto/home/users/yoichiro/software/miniconda3/envs/20220601_CB_AM_PHD2/lib/R/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
