
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
│   └── sinto.zip : customized version of sinto (installed via pip)
│       └── Modified files: fragments.py, utils.py, cli.py, arguments.py
│           Edits (denoted #AMANDA) enabled direct QNAME output during fragment creation.
│           The modified package was originally located in: /miniforge/envs/downsampling/lib/python3.10/site-packages/sinto/
│ 
└─── 02_arrows/
│   │   02_arrows.R : R script to generate 24 arrow files (via ArchR)
│   │   ArchR-createArrows*.log : logs from arrow file generation
│   │   QualityControl.tar.gz : per-replicate QC images from ArchR 
│
│
└─── 03_archr/
│   │   03_archr.R : creates ArchR project, process cell lines together, and performs quality filtering
│   │
│   └─── ENCODE_snATAC.tar.gz : ArchR project (not included due to size; downloadable separately, ~360MB)
│
│
└─── 03_filtered_bams/
│   │   03_filtered_bams.sh : filters BAMs to remove (1) low-quality reads filtered by sinto prior to fragment creation (via samtools)
|   |                                                (2) reads not mapping to chr1-22, X, or Y (via samtools)
|   |                                                (3) reads associated with filtered-out cell barcodes in ArchR (via sinto)
|   | 
│   │   *.bam : (not included due to size; downloadable separately, 568G)
|
│
└─── 03_peakcalls/ 
│   │   MACS_*.R : Rscript to (1) call peaks from BAM file using macs3 (2) standardize regions to 500bp (3) perform iterative peak overlap (adotpted from ArchR)
│   │              **NOTE** a user-friendly version is available at https://github.com/aseveritt/scBAMpler
│   │   *_500bp.exclusion.bed.gz : standardized peak sets per cell line 
│   │   Union_filt_500bp.exclusion.bed.gz : union of standardized peaks across cell lines, also processed via iterative overlap;
│   │                                       used in the cell-homogeneity and differential binding module analyses  
│   │ 
│   └─── qsubs/ : stdin/stdout logs and metadata from peak-calling job submissions
|
|
└─── 04_scripts/
│   │   scBAMpler_wrapper.ipynb : Script used to generate all scBAMpler commands and join summary statistics into one file. 
│   │
│   └───
|
|
└─── 04_downsampling/
│   └─── 00_celldicts/ : Location where all scBAMpler cell line dictionaries are stored. 
│   └─── 00_cellsim/   : Location input files for the scBAMpler cell-similarity extension are stored. 
│   └─── 01_original/  : Fragment files for 5 cell lines of interest (required for PRINT)
│   └─── 02_cells/     : Location where all cell-downsampled files are stored.
│   │   |                **Note**, the BAM and fragment files are no longer stored.
│   │   |                Only the read-QNAMEs are, which can be used to regenerate BAMs using scBAmpler/samtools 
│   │   └─── c1000/
│   │   └─── c2500/
│   │   └─── c5000/
│   │   └─── c7500/
│   │   └─── c10000/
│   │   └─── c20000/
│   │   └─── c30000/
│   │   └─── c40000/
│   │   └─── c50000/
│   │ 
│   └─── 03_reads/
│   │ 
│   └─── 04_frip/
│   │ 
│   └─── 06_cellsim/
│   │ 
│   │   sampling_stats.csv : Summary information for all the data-quality downsampled datasets. 
│   │   sampling_stats_cellsim.csv : Summary information for all the cell-homogeneity downsampled datasets. 
|
|
└─── 05_scripts/
```


---------------

