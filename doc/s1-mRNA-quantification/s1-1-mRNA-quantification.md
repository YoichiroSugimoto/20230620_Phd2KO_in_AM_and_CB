---
title: "s1-1 mRNA quantification"
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

The abundance of mRNA will be quantified by `Salmon`.


# Setup



```r
project.dir <- "/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB"

renv::restore("/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/R")
```

```
## * The library is already synchronized with the lockfile.
```

```r
processors <- 7

temp <- sapply(list.files(
    file.path(project.dir, "R/functions"),
    full.names = TRUE
), source)
```

## The input and output files and the directories


```r
annot.dir <- file.path(project.dir, "annotation/")
annot.mm39.dir <- file.path(annot.dir, "mm39_annotation")
annot.mm39.raw.dir <- file.path(annot.mm39.dir, "raw")
salmon.index.dir <- file.path(annot.dir, "mm39_annotation/salmon_indices/") 

fq.dir <- file.path(project.dir, "data/fastq")

results.dir <- file.path(project.dir, "results")
s1.dir <- file.path(results.dir, "s1")

create.dirs(c(
    results.dir,
    s1.dir
))
```

# Quantify mRNAs



```r
sample.dt <- file.path(
    project.dir, "data/sample_data/20210401_sample_data.csv"
) %>% fread

quantmRNAs <- function(dt.idx, sample.dt, s1.dir){
    
    fq.file.name <- sample.dt[dt.idx, file_name]
    sample.name <- sample.dt[dt.idx, sample_name]

    fqs <- file.path(
        fq.dir,
        paste0(fq.file.name, c("_1.fastq.gz", "_2.fastq.gz"))
    )

    salmon.dir <- file.path("/fast/AG_Sugimoto/home/users/yoichiro/software/from_source/salmon-1.9.0_linux_x86_64")
    
    paste(
        file.path(salmon.dir, "bin", "salmon"),
        "quant",
        "-p", processors,
        "-i", file.path(salmon.index.dir, "transcripts_index"),
        "-l", "IU",
        "-1", fqs[1],
        "-2", fqs[2],
        "--validateMappings",
        "--seqBias", "--gcBias", "--posBias",
        "-o", file.path(s1.dir, sample.name)
    ) %>%
        system.cat

    return()
}

lapply(
    1:nrow(sample.dt),
    quantmRNAs,
    sample.dt = sample.dt,
    s1.dir = s1.dir
)
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
## 
## [[10]]
## NULL
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
##  package     * version date (UTC) lib source
##  BiocManager   1.30.18 2022-05-18 [1] CRAN (R 4.2.1)
##  bslib         0.4.0   2022-07-16 [1] CRAN (R 4.2.1)
##  cachem        1.0.6   2021-08-19 [1] CRAN (R 4.2.1)
##  cli           3.4.1   2022-09-23 [1] CRAN (R 4.2.1)
##  colorspace    2.0-3   2022-02-21 [1] CRAN (R 4.2.1)
##  data.table  * 1.14.4  2022-10-17 [1] CRAN (R 4.2.1)
##  digest        0.6.30  2022-10-18 [1] CRAN (R 4.2.1)
##  dplyr       * 1.0.10  2022-09-01 [1] CRAN (R 4.2.1)
##  evaluate      0.17    2022-10-07 [1] CRAN (R 4.2.1)
##  fansi         1.0.3   2022-03-24 [1] CRAN (R 4.2.1)
##  fastmap       1.1.0   2021-01-25 [1] CRAN (R 4.2.1)
##  generics      0.1.3   2022-07-05 [1] CRAN (R 4.2.1)
##  ggplot2     * 3.3.6   2022-05-03 [1] CRAN (R 4.2.1)
##  glue          1.6.2   2022-02-24 [1] CRAN (R 4.2.1)
##  gtable        0.3.1   2022-09-01 [1] CRAN (R 4.2.1)
##  htmltools     0.5.3   2022-07-18 [1] CRAN (R 4.2.1)
##  jquerylib     0.1.4   2021-04-26 [1] CRAN (R 4.2.1)
##  jsonlite      1.8.2   2022-10-02 [1] CRAN (R 4.2.1)
##  khroma      * 1.9.0   2022-06-18 [1] CRAN (R 4.2.1)
##  knitr       * 1.40    2022-08-24 [1] CRAN (R 4.2.1)
##  lifecycle     1.0.3   2022-10-07 [1] CRAN (R 4.2.1)
##  magrittr    * 2.0.3   2022-03-30 [1] CRAN (R 4.2.1)
##  munsell       0.5.0   2018-06-12 [1] CRAN (R 4.2.1)
##  pillar        1.8.1   2022-08-19 [1] CRAN (R 4.2.1)
##  pkgconfig     2.0.3   2019-09-22 [1] CRAN (R 4.2.1)
##  R6            2.5.1   2021-08-19 [1] CRAN (R 4.2.1)
##  renv          0.16.0  2022-09-29 [1] CRAN (R 4.2.1)
##  rlang         1.0.6   2022-09-24 [1] CRAN (R 4.2.1)
##  rmarkdown   * 2.17    2022-10-07 [1] CRAN (R 4.2.1)
##  sass          0.4.2   2022-07-16 [1] CRAN (R 4.2.1)
##  scales        1.2.1   2022-08-20 [1] CRAN (R 4.2.1)
##  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.1)
##  stringi       1.7.8   2022-07-11 [1] CRAN (R 4.2.1)
##  stringr     * 1.4.1   2022-08-20 [1] CRAN (R 4.2.1)
##  tibble        3.1.8   2022-07-22 [1] CRAN (R 4.2.1)
##  tidyselect    1.2.0   2022-10-10 [1] CRAN (R 4.2.1)
##  utf8          1.2.2   2021-07-24 [1] CRAN (R 4.2.1)
##  vctrs         0.4.2   2022-09-29 [1] CRAN (R 4.2.1)
##  withr         2.5.0   2022-03-03 [1] CRAN (R 4.2.1)
##  xfun          0.34    2022-10-18 [1] CRAN (R 4.2.1)
##  yaml          2.3.6   2022-10-18 [1] CRAN (R 4.2.1)
## 
##  [1] /fast/AG_Sugimoto/home/users/yoichiro/software/miniconda3/envs/20220601_CB_AM_PHD2/lib/R/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
