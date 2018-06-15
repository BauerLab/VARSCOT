################################################################################
################################################################################
## Validate biochemical with data from SITE-seq with complete pipelines from our tool and CRISPOR

library(RColorBrewer)
library(VennDiagram)
library(ROCR)
library(ggplot2)
library(ggpubr)
library(reshape)

################################################################################
# Data
################################################################################

varscotDataPathMit <- "pipeline-comparison/result-siteseq-biochemical-mit.txt"
varscotDataPathRF <- "pipeline-comparison/result-siteseq-biochemical-rf.txt"
crisporDataPath <- "pipeline-comparison/crispor-siteseq-offtargets.txt"

################################################################################
# Search comparison
################################################################################

# Biochemical data
load(file="data-objects/offtargetBiochemicalData.RData")
load(file="data-objects/ontargetBiochemicalData.RData")

# VARSCOT
# mitData <- read.table(varscotDataPathMit, stringsAsFactors=F, header=F, sep="\t")
rfData <- read.table(varscotDataPathRF, stringsAsFactors=F, header=F, sep="\t")
# colnames(mitData) <- c("Chr", "Start", "End", "Targetsite", "Score", "Strand", "Sequence", "Mismatch_Number", "Mismatch_Positions")
colnames(rfData) <- c("Chr", "Start", "End", "Targetsite", "Score", "Strand", "Sequence", "Mismatch_Number", "Mismatch_Positions")

# CRISPOR
crisporData <- read.table(crisporDataPath, stringsAsFactors=F, header=T, sep="\t")

testDataForOverlap <- offtargetBiochemicalData[,c("Chr", "Start", "Offtarget_Sequence")]
testDataForOverlap[,2] <- as.character(testDataForOverlap[,2])
testOverlapVector <- apply(testDataForOverlap, 1, paste, collapse="_")

rfData.1 = subset(rfData, Score > 0.1)[,c("Chr","Start","Sequence")]
rfData.2 = subset(rfData, Score > 0.2)[,c("Chr","Start","Sequence")]
rfData.3 = subset(rfData, Score > 0.3)[,c("Chr","Start","Sequence")]
rfData.4 = subset(rfData, Score > 0.4)[,c("Chr","Start","Sequence")]
rfData.5 = subset(rfData, Score > 0.5)[,c("Chr","Start","Sequence")]
rfData.6 = subset(rfData, Score > 0.6)[,c("Chr","Start","Sequence")]
rfData.7 = subset(rfData, Score > 0.7)[,c("Chr","Start","Sequence")]
rfData.8 = subset(rfData, Score > 0.8)[,c("Chr","Start","Sequence")]
rfData.9 = subset(rfData, Score > 0.9)[,c("Chr","Start","Sequence")]

rfData.0[,'Start'] <- as.character(as.numeric(rfData.0[,'Start']) + 1)
rfData.1[,'Start'] <- as.character(as.numeric(rfData.1[,'Start']) + 1)
rfData.2[,'Start'] <- as.character(as.numeric(rfData.2[,'Start']) + 1)
rfData.3[,'Start'] <- as.character(as.numeric(rfData.3[,'Start']) + 1)
rfData.4[,'Start'] <- as.character(as.numeric(rfData.4[,'Start']) + 1)
rfData.5[,'Start'] <- as.character(as.numeric(rfData.5[,'Start']) + 1)
rfData.6[,'Start'] <- as.character(as.numeric(rfData.6[,'Start']) + 1)
rfData.7[,'Start'] <- as.character(as.numeric(rfData.7[,'Start']) + 1)
rfData.8[,'Start'] <- as.character(as.numeric(rfData.8[,'Start']) + 1)
rfData.9[,'Start'] <- as.character(as.numeric(rfData.9[,'Start']) + 1)

rfData.1 <- apply(rfData.1, 1, paste, collapse="_")
rfData.2 <- apply(rfData.2, 1, paste, collapse="_")
rfData.3 <- apply(rfData.3, 1, paste, collapse="_")
rfData.4 <- apply(rfData.4, 1, paste, collapse="_")
rfData.5 <- apply(rfData.5, 1, paste, collapse="_")
rfData.6 <- apply(rfData.6, 1, paste, collapse="_")
rfData.7 <- apply(rfData.7, 1, paste, collapse="_")
rfData.8 <- apply(rfData.8, 1, paste, collapse="_")
rfData.9 <- apply(rfData.9, 1, paste, collapse="_")

overlap.1 <- list(testOverlapVector, rfData.1)
overlap.2 <- list(testOverlapVector, rfData.2)
overlap.3 <- list(testOverlapVector, rfData.3)
overlap.4 <- list(testOverlapVector, rfData.4)
overlap.5 <- list(testOverlapVector, rfData.5)
overlap.6 <- list(testOverlapVector, rfData.6)
overlap.7 <- list(testOverlapVector, rfData.7)
overlap.8 <- list(testOverlapVector, rfData.8)
overlap.9 <- list(testOverlapVector, rfData.9)

names(overlap.1) <- c("SITE-seq", "0.1")
names(overlap.2) <- c("SITE-seq", "0.2")
names(overlap.3) <- c("SITE-seq", "0.3")
names(overlap.4) <- c("SITE-seq", "0.4")
names(overlap.5) <- c("SITE-seq", "0.5")
names(overlap.6) <- c("SITE-seq", "0.6")
names(overlap.7) <- c("SITE-seq", "0.7")
names(overlap.8) <- c("SITE-seq", "0.8")
names(overlap.9) <- c("SITE-seq", "0.9")

venn.diagram(overlap.1, "~/Desktop/0.1.png", imagetype="png")
venn.diagram(overlap.2, "~/Desktop/0.2.png", imagetype="png")
venn.diagram(overlap.3, "~/Desktop/0.3.png", imagetype="png")
venn.diagram(overlap.4, "~/Desktop/0.4.png", imagetype="png")
venn.diagram(overlap.5, "~/Desktop/0.5.png", imagetype="png")
venn.diagram(overlap.6, "~/Desktop/0.6.png", imagetype="png")
venn.diagram(overlap.7, "~/Desktop/0.7.png", imagetype="png")
venn.diagram(overlap.8, "~/Desktop/0.8.png", imagetype="png")
venn.diagram(overlap.9, "~/Desktop/0.9.png", imagetype="png")

overlapList <- list(rfOverlapVector, crisporOverlapVector, testOverlapVector, elevationOverlapVector)
names(overlapList) <- c("VARSCOT", "CRISPOR", "SITE-seq", "Elevation")

venn.diagram(overlapList, "images/siteseqSearchComparison.png", fill=brewer.pal(9, "Blues")[c(1,3,6,9)],
             imagetype="png", main="", fontfamily = 3, main.fontfamily = 3, cat.fontfamily = 3, euler.d=FALSE, scaled = FALSE)


# Find general overlap between validation data and tools
rfDataForOverlap <- rfData[,c("Chr", "Start", "Sequence")]
rfDataForOverlap[,"Start"] <-  as.character(as.numeric(rfDataForOverlap[,"Start"]) + 1)
rfOverlapVector <- apply(rfDataForOverlap, 1, paste, collapse="_")

crisporDataForOverlap <- crisporData[,c("chrom", "start", "offtargetSeq")]
crisporDataForOverlap[,2] <- as.character(crisporDataForOverlap[,2])
crisporOverlapVector <- apply(crisporDataForOverlap, 1, paste, collapse="_")


# Metrics for paper
# Missed off-targets from mismatches
missedOfftargets <- setdiff(testOverlapVector, rfOverlapVector)
indexMissedOfftargets <- sapply(missedOfftargets, function(x) which(testOverlapVector == x))

mismatchNumberTooHigh <- length(which(offtargetBiochemicalData[indexMissedOfftargets,"NM"] > 8))

# # Missed off-targets from PAM
# otherPAM <- length(missedOfftargets) - mismatchNumberTooHigh
# 
# # RF probability of additionally found off-targets
# additionalOfftargets <- setdiff(rfOverlapVector, testOverlapVector)
# 
# indexAdditionalOfftargets <- sapply(additionalOfftargets, function(x) which(rfOverlapVector == x))
# # save(indexAdditionalOfftargets, file="data-objects/indexAdditionalOfftargets.RData")
# 
# load(file = "data-objects/indexAdditionalOfftargets.RData")
# indexAlmostActive <- length(which(rfData[unlist(indexAdditionalOfftargets),"Score"] >= quantile(rfDataForBoxplot[which(rfDataForBoxplot[,2]==1),1],0.15) & rfData[unlist(indexAdditionalOfftargets),"Score"] < 0.5))
# indexAlmostActive/length(additionalOfftargets)

# Elevation
elevationFiles <- list.files("siteseq-data/elevation-search/", "\\.txt", full.names=T)
elevationDataForOverlap <- c()
for (f in elevationFiles)
{
    data <- read.table(f, header=F, stringsAsFactors=F, skip=28, sep=";")
    for (i in 1:nrow(data))
    {
        splitRecord <- strsplit(data[i,], "\t")[[1]]
        if (length(splitRecord) > 1)
        {
            if (splitRecord[7] == "+")
            {
                start <- as.numeric(splitRecord[3])
            }
            else
            {
                start <- as.numeric(splitRecord[3])
            }
            elevationDataForOverlap <- rbind(elevationDataForOverlap, c(paste("chr", splitRecord[2], sep=""), start, gsub("_","",strsplit(splitRecord[5], "/")[[1]][2])))
        }
    }
}
elevationOverlapVector <- apply(elevationDataForOverlap, 1, paste, collapse="_")


overlapList <- list(rfOverlapVector, crisporOverlapVector, testOverlapVector, elevationOverlapVector)
names(overlapList) <- c("VARSCOT", "CRISPOR", "SITE-seq", "Elevation")

venn.diagram(overlapList, "images/siteseqSearchComparison.png", fill=brewer.pal(9, "Blues")[c(1,3,6,9)],
            imagetype="png", main="", fontfamily = 3, main.fontfamily = 3, cat.fontfamily = 3, euler.d=FALSE, scaled = FALSE)

venn.diagram(overlapList, "../doc/thesis/images/siteseqSearchComparison.png", fill=brewer.pal(9, "Blues")[c(1,3,6,9)],
            imagetype="png", main="", fontfamily = 3, main.fontfamily = 3, cat.fontfamily = 3, euler.d=FALSE, scaled = FALSE)

offtargetsFound <- intersect(rfOverlapVector, testOverlapVector)

offtargetsFoundInfo <- lapply(offtargetsFound, function(x) strsplit(x, "_")[[1]])

indexFound <- c()
for (i in 1:length(offtargetsFoundInfo))
{
    for (j in 1:nrow(offtargetBiochemicalData))
    {
        if (offtargetsFoundInfo[[i]][1] == offtargetBiochemicalData[j,"Chr"] && as.numeric(offtargetsFoundInfo[[i]][2]) == offtargetBiochemicalData[j,"Start"] && offtargetsFoundInfo[[i]][3] == offtargetBiochemicalData[j,"Offtarget_Sequence"])
        {
            indexFound <- c(indexFound, j)
            break
        }
    }
}

scoresFound <- offtargetBiochemicalData[indexFound,"Score"]
scoresNotFound <- offtargetBiochemicalData[-indexFound,"Score"]

allScores <- list(scoresFound, scoresNotFound)
names(allScores) <- c("Found","Not found")

scoresForPlot <- melt(allScores)
colnames(scoresForPlot) <- c("Score", "Offtargets")

colHist <- c("springgreen4", "dodgerblue3")

png("images/histogramOfftargetsSiteSeq.png", width=700,height=400)
ggplot(scoresForPlot, aes(Score, fill=Offtargets, col=Offtargets)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity', binwidth=1) + scale_color_manual(values=colHist)+
  scale_fill_manual(values=colHist) + theme(legend.text=element_text(size=16)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

# png("../doc/thesis/images/histogramOfftargetsSiteSeq.png", width=1000,height=600)
# ggplot(scoresForPlot, aes(Score, fill=Offtargets, col=Offtargets)) + geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity', binwidth=1) + scale_color_manual(values=colHist)+
#   scale_fill_manual(values=colHist) + theme(legend.text=element_text(size=16)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
# dev.off()
#
