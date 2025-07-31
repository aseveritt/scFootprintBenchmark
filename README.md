
# Comparative evaluation of genomic footprinting algorithms for predicting transcription factor binding sites in single-cell data.

### Abstract
Transcription factors (TFs) have millions of potential binding sites across the human genome, but only a fraction are bound in a given context. Genomic footprinting aims to identify context-specific binding sites by detecting patterns in open chromatin data. While powerful, these approaches face technical challenges, especially in single-cell applications. We developed a benchmarking framework for cell-type specific footprinting and used it to evaluate the consistency, reproducibility, and equivalency of three leading methods across data quality scenarios and as a function of cell-type similarity. Peak-level read coverage emerged as the strongest predictor of stable footprints. Motivated by limited reproducibility across tools, we built an ensemble model that improved concordance with ChIP-seq. To encourage broader adoption and development of footprinting, we provide practical guidelines for robust genomic footprinting in single-cell datasets and a roadmap for extracting deeper insights about how gene regulatory networks vary across cell types in complex tissues. 

-- Citation TBD -- 

---------------
The structure of the repository:

| Directory             | Description                              |
| --------------------- | ---------------------------------------- |
| `00_unfiltered_bams/` | Source relevant ENCODE snATAC datasets |
| `01_sinto/`           | Create fragment files and extract retained read names |
| `02_arrows/`          | Prepare arrow files for each BAM (required for ArchR) |
| `03_archr/`           | Build ArchR project and perform quality filtering  |
| `03_filtered_bams/`   | Filter BAMs to include valid cell barcodes and high-quality reads |
| `03_peakcalls/`       | Call and standardize peak locations for all cell lines |
| `04_scripts/`         | Scripts used to generate downsampled BAMs with scBAMpler |
| `04_downsampling/`    | Downsampled BAM locations |
| `05_scripts/`         | Scripts used to perform footprinting |
| `05_footprinting/`    | Raw footprinting output files |
| `06_dataquality/`     | Analysis of effect of data quality downsampling on footprinting (consistency) |



```
scFootprintBenchmark
│   README.md
│
└─── 00_unfiltered_bams/
│   │   00_unfiltered_bams.sh : wget commands to download relevant ENCODE BAM files
│
│
└─── 01_sinto/
│   │   01_sinto.sh : creates fragment.tsv.gz and extracts read QNAMEs
│   └── sinto.zip   : customized version of sinto (installed via pip)
│       └── Modified files: fragments.py, utils.py, cli.py, arguments.py
│           Edits (denoted #AMANDA) enabled direct QNAME output during fragment creation.
│           The modified package was originally located in: /miniforge/envs/downsampling/lib/python3.10/site-packages/sinto/
│ 
└─── 02_arrows/
│   │   02_arrows.R             : R script to generate 24 arrow files (via ArchR)
│   │   ArchR-createArrows*.log : logs from arrow file generation
│   │   QualityControl.tar.gz   : per-replicate QC images from ArchR 
│
│
└─── 03_archr/
│   │   03_archr.R            : creates ArchR project, process cell lines together, and performs quality filtering
│   │
│   └─── ENCODE_snATAC.tar.gz : ArchR project (downloadable via Zenodo)
│
│
└─── 03_filtered_bams/
│   │   03_filtered_bams.sh  : filters BAMs to remove (1) low-quality reads filtered by sinto prior to fragment creation (via samtools)
|   |                                                (2) reads not mapping to chr1-22, X, or Y (via samtools)
|   |                                                (3) reads associated with filtered-out cell barcodes in ArchR (via sinto)
|   | 
│   │   *.bam                : not included due to size (568G). Available upon request.
|
│
└─── 03_peakcalls/ 
│   │   MACS_*.R                           : Rscript to (1) call peaks from BAM file using macs3
│   │                                                   (2) standardize regions to 500bp
│   │                                                   (3) perform iterative peak overlap (adotpted from ArchR)
│   │                                       **NOTE** a user-friendly version is available at https://github.com/aseveritt/scBAMpler
│   │   *_500bp.exclusion.bed.gz           : standardized peak sets per cell line 
│   │   Union_filt_500bp.exclusion.bed.gz  : union of standardized peaks across cell lines, also processed via iterative overlap;
│   │                                       used in the cell-homogeneity and differential binding module analyses  
│   └─── qsubs/                            : stdin/stdout logs and metadata from peak-calling job submissions
|
|
└─── 04_scripts/
│   │   scBAMpler_wrapper.ipynb       : Script generates the SGE-qsub scripts for scBAMpler and join summary statistics into one file.
|
|
└─── 04_downsampling/
│   │   sampling_stats.csv            : seed values used in downsampling -- **allows for identical recreation of dataset using scBAMpler**
│   │ 
│   └─── 00_celldicts/                : Location where all scBAMpler cell line dictionaries are stored. 
│   └─── 00_cellsim/                  : Location input files for the scBAMpler cell-similarity extension are stored. 
│   └─── 01_original/                 : Full-depth BAM/Fragment files for 5 cell lines of interest -- available upon request
│   └─── 02_cells/                    : Location where all cell-downsampled files are stored.
│   │   |                             **Note**, the BAM and fragment files are no longer stored.
│   │   |
│   │   |    *_qnames.txt             : read-QNAMES which allow for regeneration of BAMs/Fragments using scBAMpler
│   │   |                               not included due to size (~50Gb) -- available upon request
│   │   |                                         
│   │   |    cells_sampling_stats.csv : info about cell-downsampling dataset 
│   │   |    HEPG2_c5000.sh           : example input SGE bash command to create triplicates
│   │   |    HEPG2_c5000.sh.o802527   : example output of SGE job
│   │ 
│   └─── 03_reads/                    : Location where all read-downsampled files are stored.
│   │   |    *_qnames.txt             : read-QNAMES; not included due to size (~50Gb) -- available upon request                            
│   │   |    reads_sampling_stats.csv : info about read-downsampling dataset
│   │ 
│   └─── 04_frip/                     : Location where all frip-downsampled files are stored.
│   │   |    *_qnames.txt             : frip-QNAMES; not included due to size (~300Gb) -- available upon request                            
│   │   |    frip_sampling_stats.csv  : info about frip-downsampling dataset
│   │ 
│   └─── 05_cellsim/
│   │   |    similarity_summary_stats.csv
|
|
└─── 05_scripts/
│   │   00_GenerateFootprintSubmissions.ipynb : Generates the SGE-qsub scripts required for footprinting for each downsampled dataset.
│   │   01_TimeConstraints.ipynb              : Collects statistics from STDOUT files and plots results
│   │   makePrintPeakData.R                   : Input data required by PRINT was saved for each cell line prior to footprinting to save time.
│   │   PRINT_FT.R                            : PRINT's codebase is in R, this is the script each qsub command runs (hint & tobias are command-line)
│   │   PWM_FT.R                              : monaLisa is also an R-package, this is the script each qsub command runs. 
│
│ 
|
└─── 05_footprinting/
|   │    **NOTE** due to size, its difficult to host the raw results -- though all are available upon request and below I outline the size required for each.  Instead, here I host the results for a single sample for an example case (.tar.gz)
|   │ 
|   │    memory_stats.txt   : memory requirements of each tool-file pairing
|   │
│   └─── program_inputs/            : external input files required for each program
|   │
|   │
│   └─── 01_original/
|   │   └─── hint/     : Contains *info, *.bed, *mpbs.bed, for the summary, global, and local results respectively. 
|   │   └─── print/
|   │   └─── tobias/   : Contains *_ftscore.bw, *_corrected.bw, *_results.txt, and *_TF_overviews.txt, for the per-bp signal, tn5 corrected signal, global, and local results respectively
|   │   └─── pwm/      : Contains *.bed and *_mat.txt output from monaLisa
|   │
│   └─── 02_cells/
|   │   └─── hint/     : info_stats.tar.gz (5K), denovo_beds.tar.gz (4G), mpbs_beds.tar.gz (15G)
|   │   └─── print/    : 
|   │   └─── tobias/   : bigwigs.tar.gz (300G), bindetect_global_TF_results.tar.gz (2M), bindetect_local_TF_results.tar.gz (200G)
|   │   └─── pwm/      : PWMmatch_beds.tar.gz (600G), PWMmatch_mats.tar.gz (12G)
|   │   └─── macs3/    : Peak calls for downsampled datasets (2G)
|   │
│   └─── 03_reads/
|   │   HEPG2_r1e6_example.zip : ??
|   │   └─── hint/     : info_stats.tar.gz (4K), denovo_beds.tar.gz (2G), mpbs_beds.tar.gz (8G)
|   │   └─── print/    : 
|   │   └─── tobias/   : bigwigs.tar.gz (200G), bindetect_global_TF_results.tar.gz (2M), bindetect_local_TF_results.tar.gz (120G)
|   │   └─── pwm/      : PWMmatch_beds.tar.gz (300G), PWMmatch_mats.tar.gz (5G)
|   │   └─── macs3/    : Peak calls for downsampled datasets (1G)
|   │
│   └─── 04_frip/
|   │   └─── hint/     : info_stats.tar.gz (5K), denovo_beds.tar.gz (4G), mpbs_beds.tar.gz (15G)
|   │   └─── print/    : 
|   │   └─── tobias/   : bigwigs.tar.gz (350G), bindetect_global_TF_results.tar.gz (3M), bindetect_local_TF_results.tar.gz (180G)
|   │   └─── pwm/      : PWMmatch_beds.tar.gz (300G), PWMmatch_mats.tar.gz (15G)
|   │   └─── macs3/    : Peak calls for downsampled datasets (2G)
|   │
│   └─── 05_cellsim/
|   │   └─── hint/     : info_stats.tar.gz (3K), denovo_beds.tar.gz (5G), mpbs_beds.tar.gz (20G)
|   │   └─── print/    : 
|   │   └─── tobias/                     : bigwigs.tar.gz (500G), bindetect_global_TF_results.tar.gz (2M), bindetect_local_TF_results.tar.gz (250G)
|
|
|
└─── 06_dataquality
│   │   01_MakeMatrices.ipynb            : Standardizes footprinting results into matrices
|   │   02_CalculateMetrics.ipynb        : Calculates performance metrics across conditions
|   │   03_AnalysisImages.ipynb          : Main Figures for data-quality analysis
|   │   04_TFspecific.ipynb              : Analyses about what TFs generally footprint better
|   │   05_PeakConsistency.ipynb         : Analyses about the influence of peak coverage on performance metrics.
|   │
|   │   └─── hint/
|   │   │   mats.tar.gz                   : standardized matrices. Available at Zenodo
|   │   │   └─── metrics/                 
|   │   │   │   └─── no_threshold/         : contains metrics.tar.gz (per dataset performance metrics); *_bypeak.csv (per condiditon metrics calculated over peaks); and *_bytf.csv (per condiditon metrics calculated over TFs)
|   │   │   |   │    |    metrics.tar.gz   : per dataset performance metrics
|   │   │   |   │    |    *_bypeak.csv     : per condiditon metrics calculated over peaks
|   │   │   |   │    |    *_bytf.csv       : per condiditon metrics calculated over TFs
|   │   │   |   │
|   │   │   │   └─── no_threshold_peak100/ : metrics considering peaks with coverage > 100

                    
|   │   └─── print/
|   │   │   mats.tar.gz                   : standardized matrices. Available at Zenodo
|   │   │   └─── metrics/                 
|   │   │   │   └─── thresh_03/            : metrics using 0.3 as score threshold
|   │   │   |   │    |    metrics.tar.gz   : per dataset performance metrics
|   │   │   |   │    |    *_bypeak.csv     : per condiditon metrics calculated over peaks
|   │   │   |   │    |    *_bytf.csv       : per condiditon metrics calculated over TFs
|   │   │   |   │
|   │   │   │   └─── thresh_03_peak100/    : metrics using 0.3 as score threshold and coverage > 100

|   │   └─── tobias/
|   │   │   mats.tar.gz                   : standardized matrices. Available at Zenodo

|   │   └─── tfbs_universe/
|   │   └─── cellxpeak/ 
|
└─── 06_dataquality

-- need IC content file to proceed.

```


---------------

