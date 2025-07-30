#!/usr/bin/env Rscript

library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-p", "--peak"), type="character",  default=NULL, action="store",
              help="Path to reproducible peak set", 
              dest="peak_file")
parser <- add_option(parser, c("-o", "--outdir"), type="character",  default=NULL, action="store",
              help="outputdir", 
              dest="outdir")
parser <- add_option(parser, c("-m", "--motif"), type="character",  default=NULL, action="store",
              help="Path to jaspar motif.txt", 
              dest="motif_file")
parser <- add_option(parser,  c("--method"), default="matchPWM", action="store",
              help="monalisa method", 
              dest="method")

opt <- parse_args(parser)
if (is.null(opt$peak_file)) { stop("ERROR: Missing --peak argument.") }
if (is.null(opt$outdir)) { stop("ERROR: Missing --output argument.") }
if (is.null(opt$motif_file)) { stop("ERROR: Missing --motif argument.") }

suppressPackageStartupMessages({
    library(monaLisa)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(rtracklayer)
    library(TFBSTools)
})

gr <- rtracklayer::import(opt$peak_file)
gr <- unique(gr) #important if using narrowPeak files
names(gr) <- paste0(gr@seqnames, ":", gr@ranges)

jaspar_pfm <- TFBSTools::readJASPARMatrix(opt$motif_file, matrixClass=c("PFM"))
jaspar_pwm <- TFBSTools::toPWM(jaspar_pfm, pseudocounts=0.8)
names(jaspar_pwm) = sapply(jaspar_pwm, function(x) paste0(x@name,"_", x@ID))

seqs <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
                           
            
if (opt$method == "matchPWM"){                         
    res <- monaLisa::findMotifHits(query = jaspar_pwm,
                                   subject = seqs,
                                   min.score = 6.0,
                                   method = opt$method,
                                   BPPARAM = BiocParallel::SerialParam())
}

if (opt$method == "homer2"){                         
    res <- monaLisa::findMotifHits(query = jaspar_pwm,
                                   subject = seqs,
                                   min.score = 6.0,
                                   method = opt$method,
                                   homerfile = monaLisa::findHomer("homer2"),              
                                   BPPARAM = BiocParallel::SerialParam())
}
    

res$fullname = paste0(toupper(res$pwmname), "_", res$pwmid)

m <- table(factor(seqnames(res), levels = names(gr)),
           factor(res$fullname, levels = names(jaspar_pwm))) 

#so that i can pull pwm coordinates in addition to the peak_coordinates.
seqnames_split <- strsplit(as.character(seqnames(res)), "[:-]")
chromosomes <- sapply(seqnames_split, `[[`, 1)
base_starts <- as.numeric(sapply(seqnames_split, `[[`, 2))
final_starts <- base_starts + start(res) - 1
final_ends <- base_starts + end(res) - 1
final_coord <- paste0(trimws(chromosomes), ":", 
                      trimws(format(final_starts, scientific = FALSE)), "-", 
                      trimws(format(final_ends, scientific = FALSE))
                     )
res$pwm_coord <- final_coord

tmp = as.data.frame(res[, c("pwm_coord", "fullname", "score")])
tmp <- tmp[, c("pwm_coord", "seqnames", "fullname","score")]
colnames(tmp) = c("pwm_coord", "peak_coord", "motif","score")
                      

if (endsWith(basename(opt$peak_file), "_peaks.narrowPeak")) { 
    fname1 = gsub("_peaks.narrowPeak", "_PWMmatch.bed", basename(opt$peak_file))   
    fname2 = gsub("_peaks.narrowPeak", "_PWMmatch_mat.txt", basename(opt$peak_file))  
} else if (endsWith(basename(opt$peak_file), ".bed")) { 
    fname1 = gsub(".bed", "_PWMmatch.bed", basename(opt$peak_file))   
    fname2 = gsub(".bed", "_PWMmatch_mat.txt", basename(opt$peak_file))  
} else {
    fname1 = paste0(basename(opt$peak_file), "_PWMmatch.bed")
    fname2 = paste0(basename(opt$peak_file), "_PWMmatch_mat.txt")
}
                            
write.table(tmp, 
            file = paste0(opt$outdir, fname1), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(m, 
            file = paste0(opt$outdir, fname2), 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)                           
