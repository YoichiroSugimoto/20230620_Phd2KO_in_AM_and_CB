#!/bin/sh
#$ -pe smp 16
#$ -l h_vmem=10G
#$ -m ea

## Run with the following commands
# qsub /fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/master_scripts/sub_scripts/s0_run_individual_rmds_sbatch.sh   

source /fast/home/y/ysugimo/.bashrc
conda activate 20220601_CB_AM_PHD2

cd /fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/R/s0-annotation-prepocessing

Rscript ../run_rmd.R s0-1-preparation-of-Salmon-index.rmd
