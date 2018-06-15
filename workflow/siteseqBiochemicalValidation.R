################################################################################
################################################################################
## Validate biochemical with data from SITE-seq

library(BSgenome.Hsapiens.UCSC.hg38)
library(randomForest)
library(ROCR)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(reshape)

source("evalFunctions.R")

################################################################################
# Biochemical data:
################################################################################

biochemicalDataPath <- list.files("siteseq-data/biochemical-data", pattern="\\.bed",full.names=TRUE)

biochemicalData <- lapply(biochemicalDataPath, read.table, header=FALSE, stringsAsFactors=FALSE)
names(biochemicalData) <- sapply(biochemicalDataPath, function(x) strsplit(strsplit(x, "/")[[1]][3], "_")[[1]][1])

genome <- BSgenome.Hsapiens.UCSC.hg38

# Extract sequences
# for (i in 1:length(biochemicalData))
# {
#     colnames(biochemicalData[[i]]) <- c("Chr", "Start", "End", "Targetsite", "Score", "Strand")
#     biochemicalData[[i]][,"Targetsite"] <- rep(names(biochemicalData)[i], nrow(biochemicalData[[i]]))

#     sequence <- c()
#     for (j in 1:nrow(biochemicalData[[i]]))
#     {
#         if (biochemicalData[[i]][j,"Strand"] == "+")
#         {
#             sequence <- c(sequence, as.character(genome[[biochemicalData[[i]][j,"Chr"]]][as.numeric(biochemicalData[[i]][j,"Start"]):as.numeric(biochemicalData[[i]][j,"End"])]))
#         }
#         else
#         {
#             sequence <- c(sequence, as.character(reverseComplement(genome[[biochemicalData[[i]][j,"Chr"]]][as.numeric(biochemicalData[[i]][j,"Start"]):as.numeric(biochemicalData[[i]][j,"End"])])))
#         }
#     }
#     biochemicalData[[i]] <- cbind(biochemicalData[[i]], sequence)
#     names(biochemicalData[[i]])[ncol(biochemicalData[[i]])] <- "Offtarget_Sequence"
# }


# Split into on- and off-targets
# ontargetBiochemicalData <- do.call(rbind, lapply(biochemicalData, function(x) x[1,]))
# colnames(ontargetBiochemicalData)[ncol(ontargetBiochemicalData)] <- "Target_Sequence"
# ontargetBiochemicalData[,"Target_Sequence"] <- sapply(ontargetBiochemicalData[,"Target_Sequence"], as.character)

# offtargetBiochemicalData <- do.call(rbind, lapply(biochemicalData, function(x) x[-1,]))
# offtargetBiochemicalData[,"Offtarget_Sequence"] <- sapply(offtargetBiochemicalData[,"Offtarget_Sequence"], as.character)


# Extract extended on-target sequences with flanking regions for TUSCAN
# ontargetsTuscanSequence <- c()
# for (i in 1:nrow(ontargetBiochemicalData))
# {
#     if (ontargetBiochemicalData[i,"Strand"] == "+")
#     {
#         ontargetsTuscanSequence <- c(ontargetsTuscanSequence, genome[[ontargetBiochemicalData[i,"Chr"]]][(as.numeric(ontargetBiochemicalData[i,"Start"]) - 4):(as.numeric(ontargetBiochemicalData[i,"End"]) + 3)])
#     }
#     else
#     {
#         ontargetsTuscanSequence <- c(ontargetsTuscanSequence, reverseComplement(genome[[ontargetBiochemicalData[i,"Chr"]]][(as.numeric(ontargetBiochemicalData[i,"Start"]) - 3):(as.numeric(ontargetBiochemicalData[i,"End"]) + 4)]))
#     }
# }
# names(ontargetsTuscanSequence) <- ontargetBiochemicalData[,"Targetsite"]

# ontargetsTuscanSequence <- DNAStringSet(x=ontargetsTuscanSequence)


# Write on-targets with flanking regions into FASTA
# writeXStringSet(ontargetsTuscanSequence, "siteseq-data/siteseqOntargetsFlanking.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

# Write on-targets without flanking regions into FASTA
# ontargetSequence <- c(ontargetBiochemicalData[,"Target_Sequence"])
# names(ontargetSequence) <- ontargetBiochemicalData[,"Targetsite"]
# ontargetSequence <- DNAStringSet(x=ontargetSequence)
# writeXStringSet(ontargetSequence, "siteseq-data/siteseqOntargets.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")

# Write on-targets to BED file for biochemical with our pipeline
# ontargetsForBed <- ontargetBiochemicalData[,1:6]
# ontargetsForBed[,"Start"] <- ontargetsForBed[,"Start"] - 1
# write.table(ontargetsForBed, file="siteseq-data/siteseqOntargets.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)


# Check maximum number of mismatches
# mismatchNumber <- c()
# for (i in 1:nrow(offtargetBiochemicalData))
# {
#     mismatchNumber <- c(mismatchNumber, stringDiff(ontargetBiochemicalData[which(ontargetBiochemicalData[,"Targetsite"] == offtargetBiochemicalData[i,"Targetsite"]),"Target_Sequence"], offtargetBiochemicalData[i,"Offtarget_Sequence"]))
# }

# offtargetBiochemicalData <- cbind(offtargetBiochemicalData, mismatchNumber)
# colnames(offtargetBiochemicalData)[ncol(offtargetBiochemicalData)] <- "NM"


# Add on-target sequences and TUSCAN activity
# tuscanActivity <- read.table("siteseq-data/siteseqOntargetActivity.txt", header=TRUE, stringsAsFactors=FALSE)

# activity <- sapply(offtargetBiochemicalData[,"Targetsite"], function(x) tuscanActivity[which(tuscanActivity[,"ID"] == x),"Score"])
# ontargetSequences <- sapply(offtargetBiochemicalData[,"Targetsite"], function(x) ontargetBiochemicalData[which(ontargetBiochemicalData[,"Targetsite"] == x),"Target_Sequence"])

# offtargetBiochemicalData <- cbind(offtargetBiochemicalData, ontargetSequences, activity)
# colnames(offtargetBiochemicalData)[9:10] <- c("Target_Sequence", "Target_Activity")

# offtargetBiochemicalData[,"Target_Sequence"] <- sapply(offtargetBiochemicalData[,"Target_Sequence"], as.character)

# save(offtargetBiochemicalData, file="data-objects/offtargetBiochemicalData.RData")
# save(ontargetBiochemicalData, file="data-objects/ontargetBiochemicalData.RData")

load(file="data-objects/offtargetBiochemicalData.RData")
load(file="data-objects/ontargetBiochemicalData.RData")


################################################################################
# Evaluation of Random Forest only:
################################################################################

# Check if off-targets with bulges
indexInvalid <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    if (nchar(offtargetBiochemicalData[i,"Offtarget_Sequence"]) != 23)
    {
        indexInvalid <- c(indexInvalid, i)
    }
}

# No bulges

indexCanonicalPam <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    if (substr(offtargetBiochemicalData[i,"Offtarget_Sequence"], 22, 23) == "GG")
    {
        indexCanonicalPam <- c(indexCanonicalPam, i)
    }
}

featureMatrix <- buildFeatureMatrix(offtargetBiochemicalData, 23, numeric=TRUE, validate=TRUE)

load(file="data-objects/rfClassifier.RData")

rfPrediction <- predict(rfClassifier, featureMatrix, type="prob")[,2]
rfPredictionList <- lapply(1:7, function(x) rfPrediction[which(offtargetBiochemicalData[,"Score"] == x)])

rfDataForBoxplot <- melt(rfPredictionList)
colnames(rfDataForBoxplot) <- c("Probability", "Score")

# Boxplot concentration for all predictions
png("images/boxplotConcentrationLevelRF.png", width=700, height=400)
ggplot(rfDataForBoxplot , aes(x = Score, y = Probability, group=Score)) +
       geom_boxplot(fill = "seagreen3", colour = "black") +
       scale_y_continuous(name = "Probability") +
       scale_x_discrete(name = "Concentration score", labels = 1:7, limits = 1:7) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

png("../doc/thesis/images/boxplotConcentrationLevelRF.png", width=1000, height=700)
ggplot(rfDataForBoxplot , aes(x = Score, y = Probability, group=Score)) +
       geom_boxplot(fill = "seagreen3", colour = "black") +
       scale_y_continuous(name = "Probability") +
       scale_x_discrete(name = "Concentration score", labels = 1:7, limits = 1:7) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

# Inactive and active classes for off-targets (cutoff: 5 and higher is active)
class <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    if (offtargetBiochemicalData[i,"Score"] > 4)
    {
        class <- c(class, 1)
    }
    else
    {
        class <- c(class, 0)
    }
}

# All PAMs performance
rfPredictionROC <- prediction(rfPrediction, class)
rfPerformanceROC <- performance(rfPredictionROC, measure="tpr", x.measure="fpr")
rfPerformanceAUC <- performance(rfPredictionROC, measure="auc")

png("images/rfSiteSeqPredictionAll.png", width=500, height=500)
plot(rfPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
dev.off()

# Comparison canonical PAM and NGA PAM
rfNGGPredictionROC <- prediction(rfPrediction[indexCanonicalPam], class[indexCanonicalPam])
rfNGGPerformanceROC <- performance(rfNGGPredictionROC, measure="tpr", x.measure="fpr")
rfNGGPerformanceAUC <- performance(rfNGGPredictionROC, measure="auc")

rfOtherPredictionROC <- prediction(rfPrediction[-indexCanonicalPam], class[-indexCanonicalPam])
rfOtherPerformanceROC <- performance(rfOtherPredictionROC, measure="tpr", x.measure="fpr")
rfOtherPerformanceAUC <- performance(rfOtherPredictionROC, measure="auc")

png("images/rfSiteSeqPredictionCanonical.png", width=1000, height=1000)
plot(rfNGGPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
dev.off()

png("images/rfSiteSeqPredictionNonCanonical.png", width=1000, height=1000)
plot(rfOtherPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
dev.off()

# Images for thesis
png("../doc/thesis/images/rfSiteSeqPredictionCanonical.png", width=500, height=500)
plot(rfNGGPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
dev.off()

png("../doc/thesis/images/rfSiteSeqPredictionNonCanonical.png", width=500, height=500)
plot(rfOtherPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
dev.off()

################################################################################
# Score comparison with other tools:
################################################################################

rfPredictionROC <- prediction(rfPrediction, class)
rfPerformanceROC <- performance(rfPredictionROC, measure="tpr", x.measure="fpr")
rfPerformanceAUC <- performance(rfPredictionROC, measure="auc")

# CFD score (with PAM)
cfdScoring <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    cfdScoring <- c(cfdScoring, system(paste("python score-comparison/cfd-score-calculator.py --wt ", offtargetBiochemicalData[i,"Target_Sequence"], " --off ", offtargetBiochemicalData[i,"Offtarget_Sequence"], sep=""), intern=TRUE))
}

cfdScoring <- as.numeric(sapply(cfdScoring, function(x) strsplit(x, ": ")[[1]][2]))

cfdPredictionROC <- prediction(cfdScoring, class)
cfdPerformanceROC <- performance(cfdPredictionROC, measure="tpr", x.measure="fpr")
cfdPerformanceAUC <- performance(cfdPredictionROC, measure="auc")


# MIT score (without PAM)
# First parameter on-target, second off-target
mitScoring <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    mitScoring <- c(mitScoring, system(paste("python score-comparison/mit-score.py ", substr(offtargetBiochemicalData[i,"Target_Sequence"], 1, 20), " ", substr(offtargetBiochemicalData[i,"Offtarget_Sequence"], 1, 20), sep=""), intern=TRUE))
}
mitScoring <- as.numeric(mitScoring)

mitPredictionROC <- prediction(mitScoring, class)
mitPerformanceROC <- performance(mitPredictionROC, measure="tpr", x.measure="fpr")
mitPerformanceAUC <- performance(mitPredictionROC, measure="auc")


mitDataForBoxplot <- as.data.frame(cbind(mitScoring, offtargetBiochemicalData[,"Score"]))
colnames(mitDataForBoxplot) <- c("Probability", "Score")

elevationData <- read.table("elevation_output.txt", header = TRUE, stringsAsFactors = FALSE)

elevationPredictionROC <- prediction(elevationData[,"Elevation"], class)
elevationPerformanceROC <- performance(elevationPredictionROC, measure="tpr", x.measure="fpr")
elevationPerformanceAUC <- performance(elevationPredictionROC, measure="auc")


png("images/boxplotConcentrationLevelMIT.png", width=1000, height=700)
ggplot(mitDataForBoxplot , aes(x = Score, y = Probability, group=Score)) +
       geom_boxplot(fill = "lightblue1", colour = "black") +
       scale_y_continuous(name = "MIT score", limits = c(-0.1,1.8)) +
       scale_x_discrete(name = "Concentration score", labels = 1:7, limits = 1:7) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

png("../doc/thesis/images/boxplotConcentrationLevelMIT.png", width=1000, height=700)
ggplot(mitDataForBoxplot , aes(x = Score, y = Probability, group=Score)) +
       geom_boxplot(fill = "lightblue1", colour = "black") +
       scale_y_continuous(name = "MIT score", limits = c(-0.1,1.8)) +
       scale_x_discrete(name = "Concentration score", labels = 1:7, limits = 1:7) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

png("../doc/agta/boxplotConcentrationLevelMIT.png", width=600, height=300)
ggplot(mitDataForBoxplot , aes(x = Score, y = Probability, group=Score)) +
       geom_boxplot(fill = "lightblue1", colour = "black") +
       scale_y_continuous(name = "MIT score", limits = c(-0.1,1.8)) +
       scale_x_discrete(name = "Concentration score", labels = 1:7, limits = 1:7) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()


# CCTop score
cctopScoring <- c()
for (i in 1:nrow(offtargetBiochemicalData))
{
    cctopScoring <- c(cctopScoring, cctopScoreFunction(offtargetBiochemicalData[i,"Target_Sequence"], offtargetBiochemicalData[i,"Offtarget_Sequence"]))
}

cctopPredictionROC <- prediction(cctopScoring, abs(class - 1))
cctopPerformanceROC <- performance(cctopPredictionROC, measure="tpr", x.measure="fpr")
cctopPerformanceAUC <- performance(cctopPredictionROC, measure="auc")

# Decided to exclude CCTop

png("images/biochemicalDataModelComparison.png", width=700, height=500)
plot(rfPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
plot(cfdPerformanceROC, col="Brown3", lwd=2, add=TRUE)
plot(mitPerformanceROC, col="DeepSkyBlue3", lwd=2, add=TRUE)
plot(elevationPerformanceROC, col="orange", lwd=2, add=TRUE)
lines(seq(0,1,0.01), seq(0,1,0.01), lwd=2, col="grey")
legend("bottomright", c(paste("VARSCOT - AUC = ",round(rfPerformanceAUC@y.values[[1]],2)), paste("CFD Score - AUC = ",round(cfdPerformanceAUC@y.values[[1]],2)), paste("MIT Score - AUC = ",round(mitPerformanceAUC@y.values[[1]],2)), paste("Elevation - AUC = ",round(elevationPerformanceAUC@y.values[[1]],2)), "Random"), col=c("Green4", "Brown3", "DeepSkyBlue3", "orange", "grey"), lty=1, cex=1.4, lwd=2)
dev.off()

# Image for thesis
png("../doc/thesis/images/scoreComparisonAllPAMs.png", width=1000, height=700)
plot(rfPerformanceROC, col="Green4", cex.axis=1.4, cex.lab=1.4, lwd=2)
plot(cfdPerformanceROC, col="Brown3", lwd=2, add=TRUE)
plot(mitPerformanceROC, col="DeepSkyBlue3", lwd=2, add=TRUE)
legend("bottomright", c("Random Forest", "CFD Score", "MIT Score"), col=c("Green4", "Brown3", "DeepSkyBlue3"), lty=1, cex=1.4, lwd=2)
dev.off()
