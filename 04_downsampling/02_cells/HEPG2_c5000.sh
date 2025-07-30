#!/bin/bash

#$ -S /bin/bash
#$ -wd /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/02_cells/qsubs/
#$ -j y
#$ -l mem_free=38G
#$ -l scratch=25G
#$ -pe smp 8
#$ -l h_rt=24:00:00
set -e

module load CBI miniconda3
conda activate downsampling

start=$(date +%s)
echo 'working on cell type: HEPG2_c5000_s73'
python3 /pollard/data/projects/aseveritt/encode_snatacseq/04_scripts/perform_sampling.py \
            -i /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/00_celldicts/HEPG2.pickle \
            -o /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/02_cells/c5000/HEPG2_c5000_s73.txt \
            -e cells \
            -v 5000 \
            --seed 73 \
            --nproc ${NSLOTS:-1}  \
            -b /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams/HEPG2-cellFilt.bam

echo 'working on cell type: HEPG2_c5000_s31'
python3 /pollard/data/projects/aseveritt/encode_snatacseq/04_scripts/perform_sampling.py \
            -i /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/00_celldicts/HEPG2.pickle \
            -o /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/02_cells/c5000/HEPG2_c5000_s31.txt \
            -e cells \
            -v 5000 \
            --seed 31 \
            --nproc ${NSLOTS:-1}  \
            -b /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams/HEPG2-cellFilt.bam

echo 'working on cell type: HEPG2_c5000_s39'
python3 /pollard/data/projects/aseveritt/encode_snatacseq/04_scripts/perform_sampling.py \
            -i /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/00_celldicts/HEPG2.pickle \
            -o /pollard/data/projects/aseveritt/encode_snatacseq/04_downsampling/02_cells/c5000/HEPG2_c5000_s39.txt \
            -e cells \
            -v 5000 \
            --seed 39 \
            --nproc ${NSLOTS:-1}  \
            -b /pollard/data/projects/aseveritt/encode_snatacseq/03_filtered_bams/HEPG2-cellFilt.bam

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
echo "Elapsed Time: $(($end-$start)) seconds"
