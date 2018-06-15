########################################################################################################################
########################################################################################################################
# This notebook contains description and code of the data processing for the initial model building based on data
# obtained with GUIDE-seq by Tsai, S.Q. et al. (http://www.nature.com/nbt/journal/v33/n2/full/nbt.3117.html).
#
# The machine learning problem we want to address is the classification of potential off-targets into active and
# inactive off-targets. The GUIDE-seq dataset includes 13 on-targets and in total 430 active off-targets. For that
# purpose this script builds 10 datasets containing active off-targets from the GUIDE-seq dataset and an equal number
# of inactive off-targets randomly sampled (we assume that all potential off-targets not included in the GUIDE-seq
# dataset are inactive by default).
#
# In order to obtain all potential off-targets we mapped the on-targets to the reference genome hg19 with the all mapper
# RazerS3 (https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts505) 30% error rate
# (up to 6 mismatches).

########################################################################################################################
# Read mapping (command line):
# razers3 -fl pigeonhole -ng -tc 8 -i 70 -rr 100 -m 10000000 -o razers3_guideseq.sam hg19.fa ontargetsGUIDESeq.fasta
########################################################################################################################

library(xlsx)
library(reshape)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)
source("evalFunctions.R")


########################################################################################################################
# Data:
########################################################################################################################

dataPath <- "guideseq-data/datasetGUIDESeq.xlsx"               # GUIDE-seq data
tuscanPath <- "guideseq-data/guideseqOntargetActivity.txt"     # On-target activity from TUSCAN
mappingPath <- "guideseq-data/bidir_guideseq.sam"              # Mapping data (including inactive off-targets)


########################################################################################################################
# GUIDE-seq data processing:
########################################################################################################################

# Load the dataset
dataOriginal <- read.xlsx(dataPath, sheetName = "Sheet1", stringsAsFactors = FALSE)

# Change start positions to 1-based count (end positions are the same in this case)
dataOriginal[,"Start"] <- dataOriginal[,"Start"] + 1

# Define on-target length for our analysis
seqLength <- 23

# On-targets need to be removed from the off-target dataset.
# The on-target "RNF2" is excluded from the analysis as no corresponding active off-targets are included in the dataset.
# In total we get 9 on-targets and 348 off-targets for our analysis.
# The "Target_Sequence" column in the dataset includes the NGG PAM.
# As we want consider mismatches of the first PAM position between on-target and off-target as a feature, Ns are
# replaced by the actual on-target base.

indexOntarget <- which(dataOriginal$Mismatch.Total == 0) # No off-targets with perfect matches in this data
ontargets <- dataOriginal[indexOntarget,]
data <- dataOriginal[-indexOntarget,]

# Remove RNF2 on-target
ontargets <- ontargets[-which(ontargets$Targetsite == "RNF2"),]

# Mismatches of off-targets based on extracted sequence, matches order of data
load("data-objects/guideSeqDataPlot.RData")
guideSeqDataPlot <- as.data.frame(guideSeqDataPlot)
guideSeqDataPlot[,3] <- as.numeric(guideSeqDataPlot[,3])

# Plots of activity for mismatches and PAM
mismatchActivity <- lapply(1:8, function(x) data[which(guideSeqDataPlot[,"mismatches"] == x),"GUIDE.Seq.Reads"])
mismatchActivityBoxplot <- melt(mismatchActivity)
colnames(mismatchActivityBoxplot) <- c("Activity", "Mismatches")

png("../doc/thesis/images/boxplotGuideSeqMismatches.png", width=350, height=350)
ggplot(mismatchActivityBoxplot , aes(x = Mismatches, y = Activity, group=Mismatches)) +
       geom_boxplot(fill = "lightblue1", colour = "black") +
       scale_y_continuous(name = "Activity (GUIDE-seq reads)") +
       scale_x_discrete(name = "Mismatch number", labels = 1:8, limits = 1:8) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

sizeActivity <- as.data.frame(cbind(data[,"GUIDE.Seq.Reads"], as.vector(sapply(data$Offtarget_Sequence, nchar))))
# Change for plot only
index23 <- which(sizeActivity[,2] == 23)
sizeActivity[index23,2] <- 22
colnames(sizeActivity) <- c("Activity","Length")

png("../doc/thesis/images/boxplotGuideSeqSize.png", width=350, height=350)
ggplot(sizeActivity , aes(x = Length, group=Length)) +
       geom_bar(fill = "lightblue1", colour = "black", position='identity', aes(y = ..count..), width=0.8) +
#        scale_y_discrete(name = "Count") +
       scale_x_discrete(name = "Target and off-target length", labels = c("20 bp","21 bp","23 bp"), limits = 20:22) +
       theme(axis.text=element_text(size=16), axis.title=element_text(size=18))
dev.off()

geom_histogram(alpha = 0.5, aes(y = ..count..), position = 'identity', binwidth=1

pamActivity <- as.data.frame(cbind(as.numeric(data[,"GUIDE.Seq.Reads"]), as.vector(sapply(data$Offtarget_Sequence, function(x) paste("N", substr(x, nchar(x)-1, nchar(x)), sep="")))))
colnames(pamActivity) <- c("Activity","PAM")

allPAMs <- table(pamActivity[,2])
pamActivity[,2] <- as.numeric(c(1:length(allPAMs))[as.factor(pamActivity[,2])])
pamActivity[,1] <- as.numeric(as.character(pamActivity[,1]))

png("../doc/thesis/images/boxplotGuideSeqPAM.png", width=350, height=350)
ggplot(pamActivity , aes(x = PAM, y = Activity, group=PAM)) +
       geom_boxplot(fill = "lightblue1", colour = "black") +
       scale_y_continuous(name = "Activity (GUIDE-seq reads)") +
       scale_x_discrete(name = "PAM", labels = names(allPAMs), limits = 1:8) +
       theme(axis.text=element_text(angle=90, size=16), axis.title=element_text(size=18))
dev.off()


# For our analysis we want to use on-targets with 23 bp sequence length and corresponding off-targets with mismatches
# only (no indels). We require off-targets to have a valid NGG PAM for our analysis. Therefore, we need to remove
# deviating records from the dataset.

indexDeletion <- which(sapply(data$Offtarget_Sequence, function(x) grep("-",x)) > 0)
indexDiffSize <- which(sapply(data$Offtarget_Sequence, nchar) != seqLength)
indexPamMismatch <- which(dataOriginal$X3.bp.PAM...mismatches > 0)

data <- data[-unique(c(indexDeletion, indexDiffSize, indexPamMismatch)),]

data[,"Target_Sequence"] <- sapply(data[,"Targetsite"], function(x) ontargets[which(ontargets[,"Targetsite"] == x),"Offtarget_Sequence"])


# All reported sequences from GUIDE-seq datset are compared with reference genome hg19. Reported sequences differing
# from the genome are replaced by found genomic sequence.

genome <- BSgenome.Hsapiens.UCSC.hg19


checkOntargets <- matrix(NA, nrow=nrow(ontargets), ncol=2)
colnames(checkOntargets) <- c("Name", "Found_Target")
for (i in 1:nrow(ontargets))
{
    checkOntargets[i,"Name"] <- ontargets[i,"Name"]
    if (ontargets[i,"Strand"] == "+")
    {
        # Use Offtarget_Sequence sequence as Target_Sequence has a N at 21st position (Offtarget_Sequence contains actual base)
        checkOntargets[i,"Found_Target"] <- ontargets[i,"Offtarget_Sequence"] == as.character(genome[[ontargets[i,"X.Chromosome"]]][ontargets[i,"Start"]:ontargets[i,"End"]])
    }
    else
    {
        checkOntargets[i,"Found_Target"] <- ontargets[i,"Offtarget_Sequence"] == as.character(reverseComplement(genome[[ontargets[i,"X.Chromosome"]]][ontargets[i,"Start"]:ontargets[i,"End"]]))
    }
}

table(checkOntargets[,"Found_Target"])

# Write extended sequences with flanking regions for TUSCAN
# ontargetsTuscanSequence <- c()
# for (i in 1:nrow(ontargets))
# {
#     if (ontargets[i,"Strand"] == "+")
#     {
#         ontargetsTuscanSequence <- c(ontargetsTuscanSequence, genome[[ontargets[i,"X.Chromosome"]]][(as.numeric(ontargets[i,"Start"]) - 4):(as.numeric(ontargets[i,"End"]) + 3)])
#     }
#     else
#     {
#         ontargetsTuscanSequence <- c(ontargetsTuscanSequence, reverseComplement(genome[[ontargets[i,"X.Chromosome"]]][(as.numeric(ontargets[i,"Start"]) - 3):(as.numeric(ontargets[i,"End"]) + 4)]))
#     }
# }
# names(ontargetsTuscanSequence) <- ontargets[,"Targetsite"]

# ontargetsTuscanSequence <- DNAStringSet(x=ontargetsTuscanSequence)

# writeXStringSet(ontargetsTuscanSequence, "guideseq-data/guideseqOntargetsFlanking.fasta", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")


checkOfftargets <- matrix(NA, nrow=nrow(data), ncol=2)
colnames(checkOfftargets) <- c("Name", "Found_Offtarget")
for (i in 1:nrow(data))
{
    checkOfftargets[i,"Name"] <- data[i,"Name"]
    if (data[i,"Strand"] == "+")
    {
        checkOfftargets[i,"Found_Offtarget"] <- data[i,"Offtarget_Sequence"] == as.character(genome[[data[i,"X.Chromosome"]]][data[i,"Start"]:data[i,"End"]])
    }
    else
    {
        checkOfftargets[i,"Found_Offtarget"] <- data[i,"Offtarget_Sequence"] == as.character(reverseComplement(genome[[data[i,"X.Chromosome"]]][data[i,"Start"]:data[i,"End"]]))
    }
}

table(checkOfftargets[,"Found_Offtarget"])

# 7 off-targets don't match the genomic sequence. Those sequences are replaced by the real genomic sequences and
# number of mismatches is adapted. In the following more information about the differences is provided.

adaptOfftargets <- checkOfftargets[which(checkOfftargets[,"Found_Offtarget"] == FALSE),"Name"]

differenceOfftargets <- matrix(NA, nrow=length(adaptOfftargets), ncol=6)
rownames(differenceOfftargets) <- adaptOfftargets
colnames(differenceOfftargets) <- c("Target_Sequence", "Offtarget_Sequence", "Genome_Sequence", "Mismatch_Report", "Mismatch_Target_Genome", "Mismatch_Offtarget_Genome")

for (i in adaptOfftargets)
{
    dataRecord <- data[which(data[,"Name"] == i),]

    differenceOfftargets[i,"Target_Sequence"] <- dataRecord[,"Target_Sequence"]
    differenceOfftargets[i,"Offtarget_Sequence"] <- dataRecord[,"Offtarget_Sequence"]
    if (dataRecord[,"Strand"] == "+")
    {
        differenceOfftargets[i,"Genome_Sequence"] <- as.character(genome[[dataRecord[,"X.Chromosome"]]][dataRecord[,"Start"]:dataRecord[,"End"]])
    }
    else
    {
        differenceOfftargets[i,"Genome_Sequence"] <- as.character(reverseComplement(genome[[dataRecord[,"X.Chromosome"]]][dataRecord[,"Start"]:dataRecord[,"End"]]))
    }

    differenceOfftargets[i,"Mismatch_Report"] <- dataRecord[,"Mismatch.Total"]
    differenceOfftargets[i,"Mismatch_Target_Genome"] <- stringDiff(differenceOfftargets[i,"Target_Sequence"], differenceOfftargets[i,"Genome_Sequence"])
    differenceOfftargets[i,"Mismatch_Offtarget_Genome"] <- stringDiff(differenceOfftargets[i,"Offtarget_Sequence"], differenceOfftargets[i,"Genome_Sequence"])
}

differenceOfftargets <- as.data.frame(differenceOfftargets, stringsAsFactors=FALSE)


# Off-targets are corrected and information for final datasets is extracted.

dataForSampling <- data[,c("Targetsite", "X.Chromosome", "Start", "Mismatch.Total", "Offtarget_Sequence", "Target_Sequence")]
colnames(dataForSampling) <- c("Targetsite", "Chr", "Start", "NM", "Offtarget_Sequence", "Target_Sequence")

for (i in 1:nrow(differenceOfftargets))
{
    indexRecord <- which(data[,"Name"] == rownames(differenceOfftargets)[i]) # data and dataForSampling have same record ordering
    dataForSampling[indexRecord,"Offtarget_Sequence"] <- differenceOfftargets[i,"Genome_Sequence"]
    dataForSampling[indexRecord,"NM"] <- differenceOfftargets[i,"Mismatch_Target_Genome"]
}

# A different first PAM base bet ween on-target and off-target should be treated as mismatch for our analysis. Therefore,
# we need to adapt the mismatch number in the dataset. Regression scores for on-target activity from TUSCAN are also
# integrated into the data.

for (i in 1:nrow(dataForSampling))
{
    dataForSampling[i,"NM"] <- stringDiff(dataForSampling[i,"Target_Sequence"], dataForSampling[i,"Offtarget_Sequence"])
}

# Integrate the TUSCAN on-target activity into the dataset

ontargetActivity <- read.table(tuscanPath, stringsAsFactors=FALSE, header=TRUE, row.names=1)
ontargetActivityPerOfftarget <- sapply(dataForSampling[,"Targetsite"], function(x) ontargetActivity[x,"Score"])
dataForSampling <- cbind(dataForSampling, ontargetActivityPerOfftarget)
colnames(dataForSampling)[ncol(dataForSampling)] <- "Target_Activity"

# Add class information. As GUIDE-seq off-targets are all active they all get the label 1
dataForSampling <- cbind(dataForSampling, rep(1, nrow(dataForSampling)))
colnames(dataForSampling)[ncol(dataForSampling)] <- "Class"

########################################################################################################################
# Mapping data processing:
########################################################################################################################

# Load the mapping data from RazerS3
mappingData <- read.table(mappingPath, stringsAsFactors=FALSE, header=FALSE)
colnames(mappingData) <- c("Targetsite","Flag","Chr","Start","Mapq","Cigar","Rnext","Pnext","Tlen","Seq","Qual","NM","MD")

mappingData <- unique(mappingData) # Remove any potential doubled records


# Find all off-targets and on-targets from GUIDESeq data in the mapping data and delete them as they are not considered
# inactive.

# indexGuideSeq <- rep(-1, nrow(data))
# for (i in 1:nrow(data))
# {
#     for (j in 1:nrow(mappingData))
#    {
#        if (data[i,"Targetsite"] == mappingData[j,"Targetsite"] && data[i,"X.Chromosome"] == mappingData[j,"Chr"] && data[i,"Start"] == mappingData[j,"Start"])
#        {
#            indexGuideSeq[i] <- j
#            break
#        }
#    }
#}
#save(indexGuideSeq, file="data-objects/indexGuideSeq.RData")
load(file="data-objects/indexGuideSeq.RData")

# Information for final datasets is extracted.

mappingDataForSampling <- mappingData[,c("Flag", "Targetsite", "Chr", "Start", "NM")]
mappingDataForSampling[,"NM"] <- sapply(mappingDataForSampling[,"NM"], function(x) strsplit(x, ":")[[1]][3])

ontargetRecordsIndex <- which(mappingDataForSampling[,"NM"]==0)
ontargetRecords <- mappingDataForSampling[ontargetRecordsIndex,]

# Verify that all 0 mismatch hits are on-targets
all(sapply(1:nrow(ontargetRecords), function(x)
{
    found <- FALSE
    for (i in 1:nrow(ontargets))
    {
        if (ontargets[i,"Targetsite"] == ontargetRecords[x,"Targetsite"]
            && ontargets[i,"X.Chromosome"] == ontargetRecords[x,"Chr"]
            && ontargets[i,"Start"] == ontargetRecords[x,"Start"])
        {
            found <- TRUE
            break
        }
    }
    return(found)
}))


# Delete on-target and active off-target records
mappingDataForSampling <- mappingDataForSampling[-c(indexGuideSeq[which(indexGuideSeq != -1)], ontargetRecordsIndex),]

# Extract sequences for mapping data, check from flag which reads are on the reverse strand.

# sequencesMappingData <- rep(NA, nrow(mappingDataForSampling))
# for (i in 1:length(sequencesMappingData))
# {
#     if (as.integer(intToBits(mappingDataForSampling[i,"Flag"]))[5] == 0)
#     {
#         sequencesMappingData[i] <- as.character(genome[[mappingDataForSampling[i,"Chr"]]][mappingDataForSampling[i,"Start"]:(mappingDataForSampling[i,"Start"]+22)])
#     }
#     else
#     {
#         sequencesMappingData[i] <- as.character(reverseComplement(genome[[mappingDataForSampling[i,"Chr"]]][mappingDataForSampling[i,"Start"]:(mappingDataForSampling[i,"Start"]+22)]))
#     }

# }
# save(sequencesMappingData, file="data-objects/sequencesMappingData.RData")
load(file="data-objects/sequencesMappingData.RData")

# Include on-target sequences and activities into data
ontargetSequencesTable <- unique(dataForSampling[,c("Targetsite", "Target_Sequence", "Target_Activity")])
ontargetSequences <- sapply(mappingDataForSampling[,"Targetsite"], function(x) ontargetSequencesTable[which(ontargetSequencesTable[,"Targetsite"] == x), "Target_Sequence"])
ontargetActivitiesMapping <- sapply(mappingDataForSampling[,"Targetsite"], function(x) ontargetSequencesTable[which(ontargetSequencesTable[,"Targetsite"] == x), "Target_Activity"])

mappingDataForSampling <- cbind(mappingDataForSampling, sequencesMappingData, ontargetSequences, ontargetActivitiesMapping)
colnames(mappingDataForSampling)[(ncol(mappingDataForSampling)-2):ncol(mappingDataForSampling)] <- c("Offtarget_Sequence", "Target_Sequence", "Target_Activity")

# Delete flag column for further analysis
mappingDataForSampling <- mappingDataForSampling[,-1]

# Add class information. As off-targets from mapping are considered inactive (we filtered out all potential active
# off-targets) they all get the label 0.
mappingDataForSampling <- cbind(mappingDataForSampling, rep(0, nrow(mappingDataForSampling)))
colnames(mappingDataForSampling)[ncol(mappingDataForSampling)] <- "Class"


########################################################################################################################
# Downsampling and dataset building
########################################################################################################################

# With more than 1 million potential inactive off-targets the number of inactive off-targets exceeds the number of active off-targets
# multiple times. Therefore, we downsample the inactive off-targets in order to get the same number of inactive and
# active off-targets for our datasets. Downsampling is performed according to the distribution of total number of off-targets
# per on-target and a weight is applied for mismatch number to favor a small mismatch number.
# We build 10 different datasets which all contain the same active off-targets but different samples of inactive
# off-targets in order to verify in the classification process that we don't hit the outcome by chance based on a lucky
# sampling run.

# Compute the total number of off-targets per on-target that need to be sampled
# totalNumberPotentialOffTargets <- table(rbind(dataForSampling, mappingDataForSampling)[,"Targetsite"])
totalNumberPotentialOffTargets <- table(dataForSampling[,"Targetsite"])

# Divide the inactive off-targets into groups accorindg to the on-target they belong to

groupedMappingDataForSampling <- lapply(names(totalNumberPotentialOffTargets), function(x) mappingDataForSampling[which(mappingDataForSampling[,"Targetsite"] == x),])
names(groupedMappingDataForSampling) <- names(totalNumberPotentialOffTargets)

# Mismatch weights for sampling
mismatchWeights <- c(100000, 50000, 1000, 500, 100, 1, 1)
names(mismatchWeights) <- sapply(2:8, as.character)

# Define the sample weights for each group of off-targets by mismatch distribution of active dataset

sampleWeights <- lapply(groupedMappingDataForSampling, function(x) mismatchWeights[x[,"NM"]])
sampleWeights <- lapply(sampleWeights, function(x) return(x/sum(x)))

set.seed(42)


datasetsSampling <- lapply(1:10, function(x) dataForSampling)
for (i in 1:length(datasetsSampling))
{
    allSamples <-c()
    for (j in names(totalNumberPotentialOffTargets))
    {
        weights <- sampleWeights[j][[1]]
        data <- groupedMappingDataForSampling[j][[1]]
        inactiveSample <- sample(1:nrow(data), totalNumberPotentialOffTargets[j], replace=FALSE, prob=weights)
        allSamples <- rbind(allSamples, data[inactiveSample,])
    }
    datasetsSampling[[i]] <- rbind(datasetsSampling[[i]], allSamples)
}

lapply(datasetsSampling, function(x) table(x[,"NM"]))

# save(datasetsSampling, file="data-objects/datasetsSampling.RData")
