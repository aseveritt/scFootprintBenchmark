#!/usr/bin/env Rscript

#TO GENERATE RETICULATE PART
#conda_create("r-reticulate", python_version = "3.9")
#conda_install(envname ="r-reticulate", packages=c("tensorflow", "scipy","h5py"))


library(reticulate)
use_condaenv("r-reticulate")

suppressPackageStartupMessages({
    library(optparse)
    library(gtools)
    library(dplyr)
    library(TFBSTools)
    library(keras)
    library(doSNOW)
    library(hdf5r)
    library(tensorflow)
    library(pbmcapply)
    library(stringr)
    library(collapse)
    library(cladoRcpp)
})
 

parser <- OptionParser()

parser <- add_option(parser, c("-f", "--frag"), type="character", default=NULL, action="store",
              help="Concatenated fragment file",
              dest="frag_file")
parser <- add_option(parser, c("-b", "--barcodegroup"), type="character", default=NULL, action="store",
              help="Path to barcode group file",
              dest="barcode_file")
parser <- add_option(parser, c("-p", "--peak"), type="character",  default=NULL, action="store",
              help="Path to reproducible peak set", 
              dest="peak_file")
parser <- add_option(parser, c("-m", "--motif"), type="character",  default=NULL, action="store",
              help="Path to jaspar motif.txt", 
              dest="motif_file")
parser <- add_option(parser,  c("-o", "--outdir"), type="character",  default=NULL, action="store",
              help="Directory where files should be saved.", 
              dest="outdir")
parser <- add_option(parser,  c("-s", "--scale"),   default=40, action="store",
              help="number of process available to ArchR", 
              dest="scale")
parser <- add_option(parser,  c("-n", "--nproc"),   default=5, action="store",
              help="number of process available to ArchR", 
              dest="nproc")

opt <- parse_args(parser)
if (is.null(opt$frag_file)) { stop("ERROR: Missing --frag argument.") }
#if (is.null(opt$barcode_file)) { stop("ERROR: Missing --barcodegroup argument.") } #now optional
if (is.null(opt$peak_file)) { stop("ERROR: Missing --peak argument.") }
if (is.null(opt$motif_file)) { stop("ERROR: Missing --motif argument.") }
if (is.null(opt$outdir)) { stop("ERROR: Missing --outdir argument.") }
#opt$outdir <- paste0(opt$outdir, "/") #just in case. 


#CREATE OUTPUT DIRECTORY
if (!file.exists(opt$outdir)) {  dir.create(opt$outdir, recursive = T); print("Creating output directory")}

print("loading code files")
print_dir = "/pollard/data/projects/aseveritt/footprint_programs/PRINT/"
suppressPackageStartupMessages({

    for (file in c("utils.R", "getCounts.R", "getBias.R", "getFootprints.R",
               "getSubstructures.R", "visualization.R", "getGroupData.R", 
               "footprintTracking.R","getTFBS.R")){
        source(paste0(print_dir, "code/", file))
    }
})

################################################################################################
#MAKE PROJECT

project <- footprintingProject(projectName = "tmp", refGenome = "hg38")
l <- readRDS(opt$peak_file)
regionRanges(project) <- l[[2]]
regionBias(project) <- l[[1]]

################################################################################################
#LOAD KERNEL SIZES

print("loading kernel sizes") #updated from version 1. 
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-  readRDS(paste0(print_dir, "data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}
TFBindingModel(project) <- loadTFBSModel(paste0(print_dir, "data/TFBSPrediction/TFBS_model.h5"))

################################################################################################
#LOAD PWMS

print("loading pwms")
jaspar_pfm <- TFBSTools::readJASPARMatrix(opt$motif_file, matrixClass=c("PFM"))
jaspar_pwm <- TFBSTools::toPWM(jaspar_pfm, pseudocounts=0.8)

################################################################################################
#LOAD BARCODE GROUPS

print("loading barcode files")
barcodeGroups = data.table::fread(cmd=paste("zcat", opt$frag_file), header=F, sep="\t", select = 4)
barcodeGroups$V2 = 1 #setting all fragments to the same group
colnames(barcodeGroups) = c("barcode", "group")
barcodeGrouping(project) <- barcodeGroups
groups(project) <- "1"
groupCellType(project) <- stringr::str_split(gsub("\\..*","", basename(opt$frag_file)), "_")[[1]][1] #works for both

#ADD OUTPUT DIR
dataDir(project) = opt$outdir

################################################################################################
#ADD COUNT TENSORS (silently)

print("starting to getCountTensor")
sink("/dev/null")
project <- getCountTensor(project, opt$frag_file, barcodeGroups, returnCombined = F, chunkSize = 2000)
sink()


################################################################################################
#DO MULTISCALE FOOTPRINTS (optional)

ft_dir = file.path(dataDir(project), "multiscale_ft/")
dir.create(ft_dir, showWarnings = FALSE)
dir.create(paste0(ft_dir, "qsubs/"), showWarnings = FALSE)

#often breaks below this. 
save(project, jaspar_pwm, file=paste0(ft_dir, "checkpoint.Rdata"))

#write a bunch of qsub files to be submitted later if needed. 
tf_list <- names(jaspar_pwm)
chunk_size = 50
inner_cores = 13
data_chunks <- split(tf_list, cut(seq_along(tf_list), chunk_size, labels = FALSE))

for (i in 1:length(data_chunks)){
  qfile = sprintf("%s/chunk_%s.sh", paste0(ft_dir, "qsubs"), i)

  l = paste(data_chunks[[i]], collapse = ",")
  cat(c("#!/bin/bash", "",
        "#$ -S /bin/bash",
        sprintf("#$ -wd %s", paste0(ft_dir, "qsubs/")),
        "#$ -j y",
        "#$ -l mem_free=10G",
        "#$ -l h_rt=24:00:00", 
        "#$ -l x86-64-v=4", 
        sprintf("#$ -pe smp %s", inner_cores),
        "",
        "set -e", "",
        "module load CBI miniconda3 r",
        sprintf("Rscript /pollard/data/projects/aseveritt/encode_snatacseq/05_scripts/process_chunk.R -o %s --cores %d --tf_list %s", ft_dir, inner_cores, l),
        "",
        "[[ -n $JOB_ID ]] && qstat -j $JOB_ID"
        ),
      sep="\n",
      file=qfile
  )
}


################################################################################################
#ADD FOOTPRINTS

print("starting to getFootprints")
sink("/dev/null")
#will give odd rowSums() error if groupInfo not correct. 
project <- getFootprints(
  project, 
  mode = as.character(opt$scale),
  nCores = opt$nproc, 
  footprintRadius = opt$scale,
  flankRadius = opt$scale)
sink()

################################################################################################
#ADD TFBS

print("starting to getTFBS")
suppressPackageStartupMessages({
    project <- getTFBS(project, motifs = jaspar_pwm, nCores = opt$nproc)
})


checkTFBS <- function(project, chunkSize=2000) {
    #so many errors occur in this section, I implemented this function to try and catch a few common ones. 
    
  TFBSDir <- file.path(dataDir(project), "chunkedTFBSResults")
  chunkIntervals <- getChunkInterval(regionRanges(project), chunkSize = chunkSize)
  starts <- chunkIntervals[["starts"]]
  chunkInds <- 1:length(starts)
  #TFBSChunkFiles <- gtools::mixedsort(list.files(TFBSDir, full.names = TRUE))
  
  contain_errors <- c()
  
  for (chunkID in chunkInds) {
    ff <- file.path(TFBSDir, paste0("chunk_", chunkID, ".rds"))
    if (file.exists(ff)) {
      TFBSChunkData <- readRDS(ff)
      
      cond1 = !is.list(TFBSChunkData)
      cond2 = any(sapply(TFBSChunkData, function(x) !is.list(x)))
      is_last_chunk <- chunkID == length(chunkInds)
      cond3 <- !is_last_chunk && length(TFBSChunkData) != chunkSize
      
      # Check if TFBSChunkData is properly formatted
      if (cond1 || cond2 || cond3) {
        contain_errors <- c(contain_errors, chunkID)
        file.remove(ff)
      }
    } else {
      contain_errors <- c(contain_errors, chunkID)
    }
  }
  return(contain_errors)
}
                         
chunkErrors = checkTFBS(project)
print(cat("Re-doing TFBS for chunks:", chunkErrors))
project <- getTFBS(project, motifs = jaspar_pwm, nCores = opt$nproc, chunkInds = chunkErrors)
                         
TFBindingSE <- getTFBindingSE(project, nCores = opt$nproc)

################################################################################################
#OUTPUT EVERYTHING

short = gsub(".frags.tsv.bgz", "", basename(opt$frag_file))
print("starting to output data")
mcols(TFBindingSE)[groupCellType(project)] <- assay(TFBindingSE) 
write.table(rowRanges(TFBindingSE), 
            paste0(opt$outdir, str_split(short, "_")[[1]][1], "_granges.bed"), 
            sep = "\t",  quote=FALSE, row.names=FALSE)


