library("ArchR")
library(ggrepel)
library(dplyr)
library(pals)

set.seed(1)
addArchRThreads(threads = 12)
addArchRGenome("hg38")

ArrowFiles <- paste0("../02_arrows/",
                c("A549_rep1.arrow", "A549_rep2.arrow", "CALU3_rep1.arrow", 
                "GM12878_rep1.arrow", "GM12878_rep2.arrow","GM12878_rep3.arrow", 
                "HEPG2_rep1.arrow", "HEPG2_rep2.arrow", "HEPG2_rep3.arrow",
                "IMR90_rep1.arrow", "K562_rep1.arrow", "K562_rep2.arrow", "K562_rep3.arrow",
                "MCF10A_rep1.arrow", "MCF7_rep1.arrow", "MCF7_rep2.arrow", "MCF7_rep3.arrow", 
                "OCILY7_rep1.arrow", "PC9_rep1.arrow",
                "SKNSH_rep1.arrow", "SKNSH_rep2.arrow", "SKNSH_rep3.arrow", "SKNSH_rep4.arrow")
                )

#Have to add doublet scores to arrows before creating the object
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, 
    knnMethod = "LSI",
    force = TRUE, 
    LSIMethod = 1
)

### ---------------------------------- ###
#create project and do dim red
### ---------------------------------- ###

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "03_archr/ENCODE_snATAC",
  copyArrows = FALSE
)

proj$barcodes <- gsub(".*#","",proj$cellNames)
proj$CellLine <- sub("_.*", "", proj$Sample)
proj

#class: ArchRProject 
#outputDirectory: 03_archr/ENCODE_snATAC 
#samples(23): A549_rep1 A549_rep2 ... SKNSH_rep3 SKNSH_rep4
#sampleColData names(1): ArrowFiles
#cellColData names(15): Sample TSSEnrichment ... barcodes CellLine
#numberOfCells(1): 309945
#medianTSS(1): 8.202
#medianFrags(1): 11406

proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = 0.1, 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 15000, 
    dimsToUse = 1:30
)

proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.1
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj <- addTSNE(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30
)

### ---------------------------------- ###
#output a bunch of plots
### ---------------------------------- ###

p1 <- plotFragmentSizes(ArchRProj = proj)
p2 <- plotTSSEnrichment(ArchRProj = proj)

p3 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "ridges")
p4 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p5 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "ridges")
p6 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p7 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP")
p8 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "DoubletEnrichment", plotAs = "ridges")

p9 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "ridges")
p10 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin")

p11 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "BlacklistRatio", plotAs = "ridges")
p12 <- plotGroups(ArchRProj = proj, groupBy = "Sample", colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin")

g1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
g2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
g3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellLine", embedding = "UMAP")

g4 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
g5 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
g6 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellLine", embedding = "TSNE")


### messy, but this creates a nicer umap/tnse ###
tmp = proj@embeddings$UMAP$df
colnames(tmp) <- c("UMAP1", "UMAP2")
tmp$col = as.vector(proj@cellColData$CellLine)
    
a = tmp %>% group_by(col) %>% summarise_all(mean)
mine_g1 = ggplot(tmp, aes(x = UMAP1, y = UMAP2, color=col)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=tol(nrow(a))) +
  geom_text_repel(data = a, aes(label = col), color="black", bg.color="white", bg.r = 0.15, size = 4) +
  theme_bw() +
  ggtitle("Prefiltering: UMAP colored by Cell Line") +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank())
 
tmp = proj@embeddings$TSNE$df
colnames(tmp) <- c("TSNE1", "TSNE2")
tmp$col = as.vector(proj@cellColData$CellLine)
    
a = tmp %>% group_by(col) %>% summarise_all(mean)
mine_g2 = ggplot(tmp, aes(x = TSNE1, y = TSNE2, color=col)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=tol(nrow(a))) +
  geom_text_repel(data = a, aes(label = col), color="black", bg.color="white", bg.r = 0.15, size = 4) +
  theme_bw() +
  ggtitle("Prefiltering: TSNE colored by Cell Line") +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank())

pdf(paste0(getOutputDirectory(proj), "/Plots/01_PreFiltering.pdf"))
p1;p2;p3
p4; p5; p6
p7; p8; p9; 
p10; p11; p12
g1; g2; g3;
g4; g5; g6;
mine_g1; mine_g2
dev.off()   

### ---------------------------------- ###
#start filtering
### ---------------------------------- ###

proj <- filterDoublets(proj, filterRatio = 1)
#Filtering 36403 cells from ArchRProject!
#	A549_rep1 : 2972 of 17241 (17.2%)
#	A549_rep2 : 1852 of 13609 (13.6%)
#	CALU3_rep1 : 140 of 3745 (3.7%)
#	GM12878_rep1 : 512 of 7158 (7.2%)
#	GM12878_rep2 : 2776 of 26180 (10.6%)
#	GM12878_rep3 : 3 of 588 (0.5%)
#	HEPG2_rep1 : 740 of 8603 (8.6%)
#	HEPG2_rep2 : 1026 of 10133 (10.1%)
#	HEPG2_rep3 : 904 of 9510 (9.5%)
#	IMR90_rep1 : 120 of 3465 (3.5%)
#	K562_rep1 : 1037 of 10185 (10.2%)
#	K562_rep2 : 1470 of 31925 (4.6%)
#	K562_rep3 : 2844 of 22784 (12.5%)
#	MCF10A_rep1 : 490 of 7006 (7%)
#	MCF7_rep1 : 286 of 5355 (5.3%)
#	MCF7_rep2 : 2822 of 16874 (16.7%)
#	MCF7_rep3 : 2996 of 20618 (14.5%)
#	OCILY7_rep1 : 1433 of 11973 (12%)
#	PC9_rep1 : 740 of 8605 (8.6%)
#	SKNSH_rep1 : 2838 of 16847 (16.8%)
#	SKNSH_rep2 : 2869 of 17012 (16.9%)
#	SKNSH_rep3 : 2764 of 19410 (14.2%)
#	SKNSH_rep4 : 2769 of 21119 (13.1%)
	
	

pt1 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C1") & proj@cellColData$CellLine == "GM12878", ])
pt2 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C2") & proj@cellColData$CellLine == "OCILY7", ])
pt3 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C3") & proj@cellColData$CellLine == "CALU3", ])
pt4 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C4") & proj@cellColData$CellLine == "PC9", ])
pt5 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C5") & proj@cellColData$CellLine == "MCF7", ])
pt6 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C6") & proj@cellColData$CellLine == "MCF10A", ])
pt7 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C7") & proj@cellColData$CellLine == "A549", ])
pt8 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C8") & proj@cellColData$CellLine == "SKNSH", ])
pt9 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C9") & proj@cellColData$CellLine == "IMR90", ])
pt10 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C10") & proj@cellColData$CellLine == "K562", ])
pt11 = rownames(proj@cellColData[proj@cellColData$Clusters %in% c("C11") & proj@cellColData$CellLine == "HEPG2", ])
total = unique(c(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11))
length(total)
#[1] 194127
#[1] 256150 #after adding SK 

proj_sub = proj[proj$cellNames %in% total, ]
#numberOfCells(1): 256150
#medianTSS(1): 8.203
#medianFrags(1): 10686

### ---------------------------------- ###
#output updated images
### ---------------------------------- ###

g1 <- plotEmbedding(ArchRProj = proj_sub, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
g2 <- plotEmbedding(ArchRProj = proj_sub, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
g3 <- plotEmbedding(ArchRProj = proj_sub, colorBy = "cellColData", name = "CellLine", embedding = "UMAP")

tmp = proj_sub@embeddings$UMAP$df
colnames(tmp) <- c("UMAP1", "UMAP2")
tmp$col = as.vector(proj_sub@cellColData$CellLine)
    
a = tmp %>% group_by(col) %>% summarise_all(mean)
mine_g1 = ggplot(tmp, aes(x = UMAP1, y = UMAP2, color=col)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=tol(nrow(a))) +
  geom_text_repel(data = a, aes(label = col), color="black", bg.color="white", bg.r = 0.15, size = 4) +
  theme_bw() +
  ggtitle("Prefiltering: UMAP colored by Cell Line") +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank())
 
tmp = proj_sub@embeddings$TSNE$df
colnames(tmp) <- c("TSNE1", "TSNE2")
tmp$col = as.vector(proj_sub@cellColData$CellLine)
    
a = tmp %>% group_by(col) %>% summarise_all(mean)
mine_g2 = ggplot(tmp, aes(x = TSNE1, y = TSNE2, color=col)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=tol(nrow(a))) +
  geom_text_repel(data = a, aes(label = col), color="black", bg.color="white", bg.r = 0.15, size = 4) +
  theme_bw() +
  ggtitle("Prefiltering: TSNE colored by Cell Line") +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text = element_blank())
 
pdf(paste0(getOutputDirectory(proj), "/Plots/02_PostFiltering.pdf"))
g1; g2; g3;
mine_g1; mine_g2
dev.off() 

#output cell barcode: cell line file
df = proj_sub@cellColData
write.csv(df, file=paste0(getOutputDirectory(proj), "/CellLineInfo/cellColData.csv"), quote=F)


### ---------------------------------- ###
#save project 
### ---------------------------------- ###

saveArchRProject(ArchRProj = proj_sub, 
    outputDirectory = "03_archr/ENCODE_snATAC_cellFilt", 
    load = FALSE,
    dropCells=TRUE)
