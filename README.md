
# Comparative evaluation of genomic footprinting algorithms for predicting transcription factor binding sites in single-cell data.

### Abstract
Transcription factors (TFs) have millions of potential binding sites across the human genome, but only a fraction are bound in a given context. Genomic footprinting aims to identify context-specific binding sites by detecting patterns in open chromatin data. While powerful, these approaches face technical challenges, especially in single-cell applications. We developed a benchmarking framework for cell-type specific footprinting and used it to evaluate the consistency, reproducibility, and equivalency of three leading methods across data quality scenarios and as a function of cell-type similarity. Peak-level read coverage emerged as the strongest predictor of stable footprints. Motivated by limited reproducibility across tools, we built an ensemble model that improved concordance with ChIP-seq. To encourage broader adoption and development of footprinting, we provide practical guidelines for robust genomic footprinting in single-cell datasets and a roadmap for extracting deeper insights about how gene regulatory networks vary across cell types in complex tissues. 

---------------
The structure of the repository:
| Directory             | Description                              |
| --------------------- | ---------------------------------------- |
| `00_unfiltered_bams/` | Source relevant ENCODE snATAC datasets |
| `01_sinto/`           | Create fragment files; extract retained read names |
| `02_arrows/`          | In preparation for ArchR, create arrow files for each BAM |
| `03_archr/`           | Store the ArchR project itself  |
| `03_filtered_bams/`   | Filter BAMs to include proper cell barcodes and reads |
| `03_peakcalls/`       | Call Peak Coordinates for all cell-line BAMs |
| `03_scripts`          | Scripts that are used in all the 03_* directories |
| `04_downsampling/`    | Call Peak Coordinates for all cell-line BAMs |

```
scFootprintBenchmark
│   README.md
│
└─── 00_unfiltered_bams/
│   │   00_unfiltered_bams.sh : wget commands to pull ENCODE links
│   
└─── 01_sinto/
│   │   01_sinto.sh : creates fragment.tsv.gz & outputs QNAMES
│   │
│   └─── sinto/ : the sinto project itself
│       │   adflkj : in-house edits to output read QNAMES directly
│ 
└─── 02_arrows/
│   │   02_arrows.R : Rscript to generate 24 arrow files
│
└─── 03_archr/
│   │   03_archr.R : creates archR project, process cell lines together, quality filtering, outputs peak matrix for cell-similarity section. 
│   │
│   └─── ENCODE_snATAC/ : ArchR project-- not directly provided due to size limitations, but can be downloaded @ (~Gb)
│
└─── 03_filtered_bams/
│   │   *.bam : not directly provided due to size limitations, but can be downloaded @ (~Gb)
|
└─── 03_peakcalls/ 
│   │   *_500bp.exclusion.bed -- standardized cell line peak files 
│   └─── logs/ -- stdout files with details
|
└─── 03_scripts/
│   │
│   └─── 
```

