
# Comparative evaluation of genomic footprinting algorithms for predicting transcription factor binding sites in single-cell data.

### Abstract
Transcription factors (TFs) have millions of potential binding sites across the human genome, but only a fraction are bound in a given context. Genomic footprinting aims to identify context-specific binding sites by detecting patterns in open chromatin data. While powerful, these approaches face technical challenges, especially in single-cell applications. We developed a benchmarking framework for cell-type specific footprinting and used it to evaluate the consistency, reproducibility, and equivalency of three leading methods across data quality scenarios and as a function of cell-type similarity. Peak-level read coverage emerged as the strongest predictor of stable footprints. Motivated by limited reproducibility across tools, we built an ensemble model that improved concordance with ChIP-seq. To encourage broader adoption and development of footprinting, we provide practical guidelines for robust genomic footprinting in single-cell datasets and a roadmap for extracting deeper insights about how gene regulatory networks vary across cell types in complex tissues. 

---------------
The structure of the repository:
.
├── 00_unfiltered_bams/
│   ├── file1.md
├── 01_sinto/
│   ├── subdirectory
│   │   └── file3.js
│   └── file4.html
├── 02_arrows/
│   ├── file1.md
├── 03_archr/
│   ├── file1.md
├── 03_filtered_bams/
│   ├── file1.md
├── 03_peakcalls/
│   ├── file1.md
├── 03_scripts
│   ├── file1.md
└── README.md

