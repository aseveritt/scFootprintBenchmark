#!/bin/bash

#$ -S /bin/bash
#$ -wd /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls//qsubs/
#$ -j y
#$ -l mem_free=10G
#$ -l scratch=25G
#$ -pe smp 10
#$ -l h_rt=24:00:00
#$ -l x86-64-v=4
set -e

module load CBI miniforge3 r
conda activate macs
Rscript /pollard/data/projects/aseveritt/encode_snatacseq/03_scripts/MACS_peaks.R \
    --bam /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams//MCF10A-cellFilt.bam \
    -o /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/ \
    --cores ${NSLOTS:-1}

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
