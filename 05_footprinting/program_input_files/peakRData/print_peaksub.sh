#!/bin/bash

#$ -S /bin/bash
#$ -wd /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/
#$ -j y
#$ -l mem_free=20G
#$ -l scratch=50G
#$ -pe smp 16
#$ -l h_rt=24:00:00
set -e

module load CBI miniforge3
module load r

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/A549_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/A549.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/CALU3_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/CALU3.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/GM12878_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/GM12878.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/HEPG2_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/HEPG2.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/IMR90_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/IMR90.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/K562_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/K562.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/MCF10A_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/MCF10A.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/MCF7_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/MCF7.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/OCILY7_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/OCILY7.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/PC9_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/PC9.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R -p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/SKNSH_filt_500bp.exclusion.bed -o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/SKNSH.printinput.RData

Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/makePrintPeakData.R \
-p /pollard/data/projects/aseveritt/encode_snatacseq/03_peakcalls/Union_filt_500bp.exclusion.bed \
-o /pollard/data/projects/aseveritt/encode_snatacseq/05_footprinting/program_input_files/peakRData/Union.printinput.RData

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"