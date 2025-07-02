
#Author: Amanda Everitt
#in R

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-b", "--bam"), type="character",  default=NULL, action="store",
                     help="Path to bamfile", 
                     dest="bam_file")
parser <- add_option(parser, c("-o", "--outdir"), type="character",  default=NULL, action="store",
                     help="Path to output diretory", 
                     dest="out_dir")
parser <- add_option(parser, c("--skipmacs"), type="logical",  action="store_true", default=FALSE,
                     help="Should macs be run or just the standardize code", 
                     dest="skipmacs")
parser <- add_option(parser, c("-s", "--summitfile"), type="character",  default=NULL, action="store",
                     help="Path to summit file if you're not running macs", 
                     dest="summit_file")
parser <- add_option(parser, c("--cores"), type="integer",  default=1, action="store",
                     help="number of cores for mclapply", 
                     dest="ncores")

opt <- parse_args(parser)

if (is.null(opt$out_dir)) { stop("ERROR: Missing --out_dir argument.") }

if (!opt$skipmacs){ if (is.null(opt$bam_file)) { stop("ERROR: Missing --bam argument.") } }

if (opt$skipmacs){ if (is.null(opt$summit_file)){ stop("ERROR: Must add summit_file argument if macs isnt being run.") }}

################################################################################################


source("03_scripts/MACS_peaks_functions.R")

if (! opt$skipmacs) { 
    if (grepl("-cellFilt.bam", opt$bam_file)){ 
      prefix = gsub("-cellFilt.bam", "", basename(opt$bam_file)) 
    } else { 
      prefix = gsub(".bam", "", basename(opt$bam_file)) 
    }
    
    call_macs(opt$bam_file, opt$out_dir, prefix) 
    standardize_summits(paste0(opt$out_dir, prefix, "_summits.bed"), opt$out_dir, peaklen=500)
} else {
    standardize_summits(opt$summit_file, opt$out_dir, peaklen=500)
}

#########################################################################

