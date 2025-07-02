
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("--outdir"), type="character",  default=NULL, action="store",
                     help="Path to output diretory", 
                     dest="outdir")
parser <- add_option(parser, c("--outfile"), type="character",  default=NULL, action="store",
                     help="Path to output diretory", 
                     dest="outfile")
parser <- add_option(parser, c("--cores"), type="integer",  default=1, action="store",
                     help="number of cores for mclapply", 
                     dest="ncores")

opt <- parse_args(parser)

source("/03_scripts/MACS_peaks_functions.R")

make_union(infiles=list.files(path = opt$outdir, pattern = "_500bp.exclusion.bed", all.files = TRUE, full.names = FALSE), 
           outdir=opt$outdir,
           outfile=opt$outfile,
           ncores=opt$ncores)

