#In R

setwd("/oak/stanford/groups/howchang/users/seankw/Multiome/")
scRNA <- readRDS("Seurat/postHarmony/eightRetinas_postUMAP")
col <- readRDS("Seurat/postHarmony/eightRetinas_cellList")

library(ArchR)
addArchRGenome("hg38")

#Get input fragment files for each sample
sample <- "LVGret1"
setwd(paste("/oak/stanford/groups/howchang/users/seankw/Multiome/ArchR/", sample, "_wd", sep="")
inputFiles <- getInputFiles(sample)[1]
names(inputFiles) <- sample
createArrowFiles(inputFiles, force = TRUE)  

#Once all samples loaded, create merged ArchR project
setwd("/oak/stanford/groups/howchang/users/seankw/Multiome/ArchR/")
ArrowFiles <- c("LVGret1_wd/LVGret1.arrow", "LVGret2_wd/LVGret2.arrow", "LGS1ODret_wd/LGS1ODret.arrow", "LGS1OSret_wd/LGS1OSret.arrow", "LGS2ODret_wd/LGS2ODret.arrow", "LGS2OSret_wd/LGS2OSret.arrow", "LGS3ODret_wd/LGS3ODret.arrow", "LGS3OSret_wd/LGS3OSret.arrow") #list all samples here
proj <- ArchRProject(ArrowFiles, copyArrows = FALSE)

#Isolate cells used in Seurat
proj <- subsetCells(ArchRProj = proj, cellNames = col)

#Integrate gene expression data
seRNA1 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LVGret1/outs/filtered_feature_bc_matrix.h5"),
    names = c("LVGret1")
)
seRNA2 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LVGret2/outs/filtered_feature_bc_matrix.h5"),
    names = c("LVGret2")
)
seRNA3 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS1ODret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS1ODret")
)
seRNA4 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS1OSret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS1OSret")
)
seRNA5 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS2ODret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS2ODret")
)
seRNA6 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS2OSret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS2OSret")
)
seRNA7 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS3ODret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS3ODret")
)
seRNA8 <- import10xFeatureMatrix(   #Repeat for each sample
    input = c("/oak/stanford/groups/howchang/users/seankw/Multiome/cellranger/LGS3OSret/outs/filtered_feature_bc_matrix.h5"),
    names = c("LGS3OSret")
)

seRNA <- cbind(seRNA1, seRNA2, seRNA3, seRNA4, seRNA5, seRNA6, seRNA7, seRNA8)  #Combine all RangedSummarizedExperiment objects

proj <- addGeneExpressionMatrix(
  input = proj,
  seRNA = seRNA,
  chromSizes = getChromSizes(proj),
  excludeChr = c("chrM", "chrY"),
  scaleTo = 10000,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addGeneExpressionMatrix")
)

#Transfer cluster and identities from scRNA
clusters <- as.character(scRNA$seurat_clusters)
identities <- as.character(scRNA@active.ident)
proj$cluster <- clusters
proj$identities <- identities

#Clean up environment
rm(scRNA, col, clusters, identities, shared, seRNA1, seRNA2, seRNA3, seRNA4, seRNA5, seRNA6, seRNA7, seRNA8)
saveArchRProject(ArchRProj = proj, outputDirectory = "Saved-Projects", load = FALSE)

#Pseudobulk ATAC and Call Peaks
test <- addGroupCoverages(ArchRProj = proj, groupBy = "identities")
pathToMacs2 <- findMacs2()
test <- addReproduciblePeakSet(
    ArchRProj = test, 
    groupBy = "identities", 
    pathToMacs2 = pathToMacs2,
    maxPeaks = 10000000
)
test <- addPeakMatrix(test)
getAvailableMatrices(test)

markersPeaks <- getMarkerFeatures(
    ArchRProj = test, 
    useMatrix = "PeakMatrix", 
    groupBy = "identities",
    useGroups = c("Rod", "OFF-cone bipolar", "Muller glia", "ON-cone bipolar", "Rod bipolar", "Cone", "GABA-amacrine", "Horizontal", "Gly-amacrine", "Retinal ganglion cell", "AII-amacrine", "Astrocyte", "Microglia"),
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

#Plot heatmap
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
plotPDF(heatmapPeaks, name = "Peak-Marker_noMax-Heatmap", width = 8, height = 6, ArchRProj = test, addDOC = FALSE)

#LSI-ATAC
test <- addIterativeLSI(
    ArchRProj = test, 
    clusterParams = list(
      resolution = 0.5,
      dimsToUse = 1:20, 
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix", 
    # depthCol = "nFrags",
    iterations = 10,
    seed = 99,
    name = "LSI_ATAC",
    force = TRUE
)

#Plot UMAP
test <- addUMAP(test, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
p1 <- plotEmbedding(test, colorBy = "cellColData", name = "identities", embedding = "UMAP_ATAC", size = 0.5, labelAsFactors=F, labelMeans=F)
plotPDF(p1, name = "UMAP", ArchRProj = test, addDOC = FALSE, width = 5, height = 5)

#Correct batch effects
test <- addHarmony(
    ArchRProj = test,
    reducedDims = "LSI_ATAC",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

#Co-accessibility analysis
test <- addCoAccessibility(
    ArchRProj = test,
    reducedDims = "Harmony"
)

#Peak2gene analysis
test <- addPeak2GeneLinks(
    ArchRProj = test,
    useMatrix = "GeneExpressionMatrix",
    reducedDims = "Harmony"
)

#Create accessibility list organized by cell type
celltypes <- c("Rod", "OFF.cone.bipolar", "Muller.glia", "ON.cone.bipolar", "Rod.bipolar", "Cone", "GABA.amacrine", "Horizontal", "Gly.amacrine", "Retinal.ganglion.cell", "AII.amacrine", "Astrocyte", "Microglia")

cellPeaksBed <- NULL
for (type in celltypes) {
    list <- readRDS(paste0("ArchROutput/PeakCalls/", type, "-reproduciblePeaks.gr.rds"))
    chr <- as.vector(list %>% seqnames(.))
    start <- list %>% start(.)
    end <- list %>% end(.)
    type <- rep(type, each=length(list))    
    typeBed <- data.frame(chr, start, end, type)
    cellPeaksBed <- rbind(cellPeaksBed, typeBed)
}
write.table(cellPeaksBed, "cellPeaks.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Create promoter co-accessibility bed file
library(dplyr)

cA <- getCoAccessibility(
    ArchRProj = test,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = FALSE
)
cA_Chr <- (metadata(cA)$peakSet %>% seqnames(.))[cA$queryHits]
cA_Start <- (metadata(cA)$peakSet %>% start(.))[cA$queryHits]
cA_End <- (metadata(cA)$peakSet %>% end(.))[cA$queryHits]
cA_SubjectStart <- (metadata(cA)$peakSet %>% start(.))[cA$subjectHits]
id <- paste(cA_Chr, cA_SubjectStart)
cA_Bed <- data.frame(cA_Chr, cA_Start, cA_End, id)

pro_Chr <- test@peakSet@seqnames    #Obtain promoter peaks
pro_Start <- test@peakSet@ranges@start
id <- paste(pro_Chr, pro_Start)
peakType <- test@peakSet$peakType
PromoterPeaks <- data.frame(id, peakType)
PromoterPeaks <- PromoterPeaks[peakType == "Promoter",]

Promoter_cA <- inner_join(cA_Bed, PromoterPeaks, by="id")   #Intersect promoter peaks with subject peaks
Promoter_cA <- Promoter_cA[, -c(5)]
write.table(Promoter_cA, "Promoter_cA.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Create peak2gene bed file
p2g <- getPeak2GeneLinks(
    ArchRProj = test,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = FALSE
)
p2gChr <- (metadata(p2g)$peakSet %>% seqnames(.))[p2g$idxATAC]
p2gStart <- (metadata(p2g)$peakSet %>% start(.))[p2g$idxATAC]
p2gEnd <- (metadata(p2g)$peakSet %>% end(.))[p2g$idxATAC]
p2gGene <- (metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2gBed <- data.frame(p2gChr, p2gStart, p2gEnd, p2gGene)
write.table(p2gBed, "p2g.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Find enriched motifs 
test <- addMotifAnnotations(ArchRProj = test, motifSet = "cisbp", name = "Motif", force = TRUE)
markerPeaks <- readRDS("ArchROutput/markerPeaks.rds")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = test,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = FALSE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap-cisbp", width = 6, height = 20, ArchRProj = test, addDOC = FALSE)
write.csv(assays(enrichMotifs)$mlog10Padj, "TF_cisbp_mlog10Padj.csv")
write.csv(assays(enrichMotifs)$Enrichment, "TF_cisbp_log2FC.csv")

#Footprinting
motifPositions <- getPositions(test)
motifs <- c("OTX2", "LHX2", "ONECUT1", "SPI1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs
seFoot <- getFootprints(
  ArchRProj = test, 
  positions = motifPositions[markerMotifs], 
  groupBy = "identities"
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = test, 
  normMethod = "Subtract", 
  plotName = "Footprints-Subtract-Bias_cisbp",
  addDOC = FALSE,
  smoothWindow = 10
)

proj <- test   
saveArchRProject(ArchRProj = proj, outputDirectory = "Saved-Projects", load = FALSE)

####################################################
#Plot browser tracks
####################################################
library(GenomicRanges)
gr=GRanges(seqnames=c("chr1"),
           ranges=IRanges(start=c(44017500-50000),
                          end=c(44017500+50000)),
           strand=c("-")
)
gr

p2g <- getPeak2GeneLinks(
    ArchRProj = test,
    corCutOff = 0.5,
    resolution = 1000,
    returnLoops = TRUE
)

p <- plotBrowserTrack(
    region = gr,  #set to gr if specifying range
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(10, 3, 1.5, 4),
    ArchRProj = test, 
    groupBy = "identities", 
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE),
    upstream = 20000,
    downstream = 20000,
    loops = p2g,
    ylim = c(0,0.99)
)
plotPDF(p, name = "SLC6A9_revision", width = 5, height = 7.5, ArchRProj = test, addDOC = FALSE)

########################################
#Export bigwig files
########################################
getGroupBW(
  ArchRProj = test,
  groupBy = "identities",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
