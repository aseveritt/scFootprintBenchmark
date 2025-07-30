suppressPackageStartupMessages({
    library(optparse)
    library(gtools) #otherwise it breaks
    library(TFBSTools)
})

parser <- OptionParser()
parser <- add_option(parser, c("-p", "--peak"), type="character",  default=NULL, action="store",
              help="Path to reproducible peak set", 
              dest="peak_file")
parser <- add_option(parser,  c("-o", "--outfile"), type="character",  default=NULL, action="store",
              help="Directory where files should be saved.", 
              dest="out_file")
opt <- parse_args(parser)


suppressPackageStartupMessages({
    print_dir = "/pollard/data/projects/aseveritt/footprint_programs/PRINT/"
    for (file in c("utils.R", "getCounts.R", "getBias.R", "getFootprints.R",
               "getSubstructures.R", "visualization.R", "getGroupData.R", 
               "footprintTracking.R","getTFBS.R")){
      source(paste0(print_dir, "code/", file))
    }
})

project <- footprintingProject(projectName = "tmp", refGenome = "hg38")
regionBed <- read.table(opt$peak_file, header = F)
regions <- GRanges(paste0(regionBed$V1, ":", regionBed$V2, "-", regionBed$V3))

#you have to normalize peaks before the step below otherwise you get a cryptic errror. 
#ended up changing the workflow prior to allow for this, leaving comment as reminder. 
#normalize_peaks = 1000
#regions <- resize(regions, normalize_peaks, fix = "center")

regionRanges(project) <- regions

setwd("/pollard/data/projects/aseveritt/footprint_programs/PRINT/analyses/TFBSPrediction")
project <- getPrecomputedBias(project, nCores = 12)

saveRDS(list(regionBias(project), regionRanges(project)), file = opt$out_file)
