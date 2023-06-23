---
title: "sp-1 Preparation of Salmon index"
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

Salmon index is prepared for the quantification of RNA-Seq data.

# Setup



```r
project.dir <- "/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB"

renv::restore("/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/R")
```

```
## * The library is already synchronized with the lockfile.
```

```r
library("rtracklayer")

options(timeout = 10000) # Downloading genome takes a long time

processors <- 15

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
annot.ps.dir <- file.path(annot.dir, "mm39_annotation/processed_data/")

create.dirs(c(
    annot.dir,
    annot.mm39.dir,
    annot.mm39.raw.dir,
    salmon.index.dir,
    annot.ps.dir
))
```

# Download annotation data



```r
gff3.url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.primary_assembly.annotation.gff3.gz"
tx.fasta.url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.transcripts.fa.gz"
genome.fasta.url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/GRCm39.primary_assembly.genome.fa.gz"

file.prefix <- "GRCm39_GENCODE_M31"

gff3.file <- file.path(
    annot.mm39.raw.dir,
    paste0(
        file.prefix,
        "_",
        format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
        ".gff3.gz")
)

download.file(
    url = gff3.url,
    destfile = gff3.file
)

tx.fasta.file <- file.path(
    annot.mm39.raw.dir,
    paste0(
        file.prefix,
        "_",
        format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
        "_tx_fasta.gz")
)

download.file(
    url = tx.fasta.url,
    destfile = tx.fasta.file
)

genome.fasta.file <- file.path(
    annot.mm39.raw.dir,
    paste0(
        file.prefix,
        "_",
        format(as.POSIXlt(Sys.time(), "GMT"), c("%Y_%m%d_%H%M%S")),
        "_genome_fasta.gz")
)

download.file(
    url = genome.fasta.url,
    destfile = genome.fasta.file
)
```

# Summarise transcript information



```r
all.tx.gr <- rtracklayer::import(gff3.file)

all.tx.dt <- as.data.frame(all.tx.gr) %>%
    data.table

all.tx.dt <- all.tx.dt[, .(
    transcript_id, gene_id, gene_name,
    gene_type, transcript_type, seqnames
)][!is.na(transcript_id)][!duplicated(transcript_id)]

all.tx.data.file <- file.path(
    annot.ps.dir,
    "all-tx-info.csv"
)

fwrite(all.tx.dt, file = all.tx.data.file)
```

# Preparation for Salmon index generation



```r
decoy.tx.file <- file.path(
    annot.ps.dir,
    "decoy.txt"
)

genome.bs <- rtracklayer::import(
                              genome.fasta.file,
                              format = "fasta"
                          )

decoy.dt <- data.table(
    seqnames = names(genome.bs)
)

decoy.dt[, seqnames := str_split_fixed(seqnames, " ", n = 2)[, 1]]

fwrite(
    decoy.dt,
    decoy.tx.file, col.names = FALSE
)

gentrome.fa.file <- file.path(
    annot.ps.dir,
    "gentrome.fa.gz"
)

temp <- paste(
    "cat",
    tx.fasta.file, genome.fasta.file, ">",
    gentrome.fa.file
) %>%
    system.cat
```


# Create Salmon index



```r
salmon.dir <- file.path("/fast/AG_Sugimoto/home/users/yoichiro/software/from_source/salmon-1.9.0_linux_x86_64")

paste(
    file.path(salmon.dir, "bin", "salmon"),
    "index",
    "-t", gentrome.fa.file,
    "-d", decoy.tx.file,
    "-p", processors,
    "-i", file.path(salmon.index.dir, "transcripts_index"),
    "--gencode"
) %>%
    system.cat
```

```
##  [1] "Threads = 15"                                                                                                                                              
##  [2] "Vertex length = 31"                                                                                                                                        
##  [3] "Hash functions = 5"                                                                                                                                        
##  [4] "Filter size = 68719476736"                                                                                                                                 
##  [5] "Capacity = 2"                                                                                                                                              
##  [6] "Files: "                                                                                                                                                   
##  [7] "/fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/annotation//mm39_annotation/salmon_indices//transcripts_index/ref_k31_fixed.fa"
##  [8] "--------------------------------------------------------------------------------"                                                                          
##  [9] "Round 0, 0:68719476736"                                                                                                                                    
## [10] "Pass\tFilling\tFiltering"                                                                                                                                  
## [11] "1\t94\t189\t"                                                                                                                                              
## [12] "2\t40\t5"                                                                                                                                                  
## [13] "True junctions count = 15001189"                                                                                                                           
## [14] "False junctions count = 1452861"                                                                                                                           
## [15] "Hash table size = 16454050"                                                                                                                                
## [16] "Candidate marks count = 339394090"                                                                                                                         
## [17] "--------------------------------------------------------------------------------"                                                                          
## [18] "Reallocating bifurcations time: 6"                                                                                                                         
## [19] "True marks count: 337697206"                                                                                                                               
## [20] "Edges construction time: 50"                                                                                                                               
## [21] "--------------------------------------------------------------------------------"                                                                          
## [22] "Distinct junctions = 15001189"                                                                                                                             
## [23] ""                                                                                                                                                          
## [24] "for info, total work write each  : 2.331    total work inram from level 3 : 4.322  total work raw : 25.000 "                                               
## [25] "Bitarray     11905327552  bits (100.00 %)   (array + ranks )"                                                                                              
## [26] "final hash        213024  bits (0.00 %) (nb in final hash 634)"
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
##  package              * version   date (UTC) lib source
##  Biobase                2.56.0    2022-04-26 [1] Bioconductor
##  BiocGenerics         * 0.42.0    2022-04-26 [1] Bioconductor
##  BiocIO                 1.6.0     2022-04-26 [1] Bioconductor
##  BiocManager            1.30.18   2022-05-18 [1] CRAN (R 4.2.1)
##  BiocParallel           1.30.4    2022-10-11 [1] Bioconductor
##  Biostrings             2.64.1    2022-08-18 [1] Bioconductor
##  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.2.1)
##  bslib                  0.4.0     2022-07-16 [1] CRAN (R 4.2.1)
##  cachem                 1.0.6     2021-08-19 [1] CRAN (R 4.2.1)
##  cli                    3.4.1     2022-09-23 [1] CRAN (R 4.2.1)
##  codetools              0.2-18    2020-11-04 [1] CRAN (R 4.2.1)
##  colorspace             2.0-3     2022-02-21 [1] CRAN (R 4.2.1)
##  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.2.1)
##  data.table           * 1.14.4    2022-10-17 [1] CRAN (R 4.2.1)
##  DelayedArray           0.22.0    2022-04-26 [1] Bioconductor
##  digest                 0.6.30    2022-10-18 [1] CRAN (R 4.2.1)
##  dplyr                * 1.0.10    2022-09-01 [1] CRAN (R 4.2.1)
##  evaluate               0.17      2022-10-07 [1] CRAN (R 4.2.1)
##  fansi                  1.0.3     2022-03-24 [1] CRAN (R 4.2.1)
##  fastmap                1.1.0     2021-01-25 [1] CRAN (R 4.2.1)
##  generics               0.1.3     2022-07-05 [1] CRAN (R 4.2.1)
##  GenomeInfoDb         * 1.32.4    2022-09-06 [1] Bioconductor
##  GenomeInfoDbData       1.2.8     2022-10-21 [1] Bioconductor
##  GenomicAlignments      1.32.1    2022-07-24 [1] Bioconductor
##  GenomicRanges        * 1.48.0    2022-04-26 [1] Bioconductor
##  ggplot2              * 3.3.6     2022-05-03 [1] CRAN (R 4.2.1)
##  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.2.1)
##  gtable                 0.3.1     2022-09-01 [1] CRAN (R 4.2.1)
##  htmltools              0.5.3     2022-07-18 [1] CRAN (R 4.2.1)
##  IRanges              * 2.30.1    2022-08-18 [1] Bioconductor
##  jquerylib              0.1.4     2021-04-26 [1] CRAN (R 4.2.1)
##  jsonlite               1.8.2     2022-10-02 [1] CRAN (R 4.2.1)
##  khroma               * 1.9.0     2022-06-18 [1] CRAN (R 4.2.1)
##  knitr                * 1.40      2022-08-24 [1] CRAN (R 4.2.1)
##  lattice                0.20-45   2021-09-22 [1] CRAN (R 4.2.1)
##  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.2.1)
##  magrittr             * 2.0.3     2022-03-30 [1] CRAN (R 4.2.1)
##  Matrix                 1.5-1     2022-09-13 [1] CRAN (R 4.2.1)
##  MatrixGenerics         1.8.1     2022-06-26 [1] Bioconductor
##  matrixStats            0.62.0    2022-04-19 [1] CRAN (R 4.2.1)
##  munsell                0.5.0     2018-06-12 [1] CRAN (R 4.2.1)
##  pillar                 1.8.1     2022-08-19 [1] CRAN (R 4.2.1)
##  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.2.1)
##  R6                     2.5.1     2021-08-19 [1] CRAN (R 4.2.1)
##  RCurl                  1.98-1.9  2022-10-03 [1] CRAN (R 4.2.1)
##  renv                   0.16.0    2022-09-29 [1] CRAN (R 4.2.1)
##  restfulr               0.0.15    2022-06-16 [1] CRAN (R 4.2.1)
##  rjson                  0.2.21    2022-01-09 [1] CRAN (R 4.2.1)
##  rlang                  1.0.6     2022-09-24 [1] CRAN (R 4.2.1)
##  rmarkdown            * 2.17      2022-10-07 [1] CRAN (R 4.2.1)
##  Rsamtools              2.12.0    2022-04-26 [1] Bioconductor
##  rtracklayer          * 1.56.1    2022-06-23 [1] Bioconductor
##  S4Vectors            * 0.34.0    2022-04-26 [1] Bioconductor
##  sass                   0.4.2     2022-07-16 [1] CRAN (R 4.2.1)
##  scales                 1.2.1     2022-08-20 [1] CRAN (R 4.2.1)
##  sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.2.1)
##  stringi                1.7.8     2022-07-11 [1] CRAN (R 4.2.1)
##  stringr              * 1.4.1     2022-08-20 [1] CRAN (R 4.2.1)
##  SummarizedExperiment   1.26.1    2022-04-29 [1] Bioconductor
##  tibble                 3.1.8     2022-07-22 [1] CRAN (R 4.2.1)
##  tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.2.1)
##  utf8                   1.2.2     2021-07-24 [1] CRAN (R 4.2.1)
##  vctrs                  0.4.2     2022-09-29 [1] CRAN (R 4.2.1)
##  withr                  2.5.0     2022-03-03 [1] CRAN (R 4.2.1)
##  xfun                   0.34      2022-10-18 [1] CRAN (R 4.2.1)
##  XML                    3.99-0.11 2022-10-03 [1] CRAN (R 4.2.1)
##  XVector                0.36.0    2022-04-26 [1] Bioconductor
##  yaml                   2.3.6     2022-10-18 [1] CRAN (R 4.2.1)
##  zlibbioc               1.42.0    2022-04-26 [1] Bioconductor
## 
##  [1] /fast/AG_Sugimoto/home/users/yoichiro/software/miniconda3/envs/20220601_CB_AM_PHD2/lib/R/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
