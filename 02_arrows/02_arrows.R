library("ArchR")
set.seed(1)
addArchRThreads(threads = 14)
addArchRGenome("hg38")

inputFiles <- paste0(
                "/pollard/data/projects/aseveritt/encode_snatacseq/01_sinto/",
                c("A549_rep1_fragments_sinto.tsv.gz", "A549_rep2_fragments_sinto.tsv.gz",
                  "CALU3_rep1_fragments_sinto.tsv.gz", "GM12878_rep1_fragments_sinto.tsv.gz", 
                  "GM12878_rep2_fragments_sinto.tsv.gz", "GM12878_rep3_fragments_sinto.tsv.gz",
                  "HEPG2_rep1_fragments_sinto.tsv.gz", "HEPG2_rep2_fragments_sinto.tsv.gz", 
                  "HEPG2_rep3_fragments_sinto.tsv.gz", "IMR90_rep1_fragments_sinto.tsv.gz",
                  "K562_rep1_fragments_sinto.tsv.gz", "K562_rep2_fragments_sinto.tsv.gz", 
                  "K562_rep3_fragments_sinto.tsv.gz","MCF10A_rep1_fragments_sinto.tsv.gz",
                  "MCF7_rep1_fragments_sinto.tsv.gz","MCF7_rep2_fragments_sinto.tsv.gz", 
                  "MCF7_rep3_fragments_sinto.tsv.gz","OCILY7_rep1_fragments_sinto.tsv.gz",
                  "PC9_rep1_fragments_sinto.tsv.gz", "SKNSH_rep1_fragments_sinto.tsv.gz", 
                  "SKNSH_rep2_fragments_sinto.tsv.gz", "SKNSH_rep3_fragments_sinto.tsv.gz",
                  "SKNSH_rep4_fragments_sinto.tsv.gz")
    )
#addArchRGenome("hg38")
mynames <- c("A549_rep1", "A549_rep2", "CALU3_rep1",
			 "GM12878_rep1", "GM12878_rep2","GM12878_rep3", 
             "HEPG2_rep1", "HEPG2_rep2", "HEPG2_rep3",
             "IMR90_rep1", 
             "K562_rep1", "K562_rep2", "K562_rep3",
             "MCF10A_rep1", 
             "MCF7_rep1", "MCF7_rep2", "MCF7_rep3", 
             "OCILY7_rep1", "PC9_rep1",
             "SKNSH_rep1", "SKNSH_rep2", "SKNSH_rep3", "SKNSH_rep4")

setwd("/pollard/data/projects/aseveritt/encode_snatacseq/02_arrows/")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  outputNames = mynames,
  sampleNames = mynames,
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  excludeChr = c("chrM")  
)
