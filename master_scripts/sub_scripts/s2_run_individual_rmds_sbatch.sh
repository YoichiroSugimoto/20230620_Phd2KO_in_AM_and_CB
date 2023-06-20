#!/bin/sh
#$ -pe smp 8
#$ -l h_vmem=10G
#$ -m ea

## Run with the following commands
# qsub /fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/master_scripts/sub_scripts/s2_run_individual_rmds_sbatch.sh   

source /fast/home/y/ysugimo/.bashrc
conda activate 20230620_Phd2KO_in_AM_and_CB

cd /fast/AG_Sugimoto/home/users/yoichiro/projects/20230620_Phd2KO_in_AM_and_CB/R/s2-differential-expression-analysis

Rscript ../run_rmd.R s2-1-differential-expression-analysis.rmd
Rscript ../run_rmd.R s2-2-master-table.rmd