#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args) != 3)
{
  stop("Error: Wrong number of arguments supplied. Arguments must be offtarget path (pipeline output), feature matrix path and probability (true/false).")
}

# Build feature matrix from pipeline output
# A numeric feature matrix is built (factors produce problems when levels are deviating between training and test data)
#
# Pipeline output: Targetsite - Chromosome - Start - Strand - Mismatch Positions - SNP Information

suppressMessages(library(randomForest))

# Function that processes the pipeline output to build a dataset
# Extracts sequences and includes on-target activity
# Must return dataset with columns: Targetsite - Chr - Start - Strand - NM - Target_Sequence - Offtarget_Sequence - Target_Activity
# Perform the classification

rfClassification <- function(offtargetPath, featureMatrixPath, prob=FALSE)
{
    load(file="classification/rfClassifier.RData")

    featureMatrix <- read.table(featureMatrixPath, stringsAsFactors=FALSE, header=TRUE, row.names=1)

    if (prob)
    {
        rfPrediction <- predict(rfClassifier, featureMatrix, type="prob")[,2]
    }
    else
    {
        rfPrediction <- predict(rfClassifier, featureMatrix)
    }

    offtargetData <- read.table(offtargetPath, stringsAsFactors=FALSE, header=FALSE, sep="\t")
    if (ncol(offtargetData) == 10)
    {
        colnames(offtargetData) <- c("#Chr", "Start", "End", "Name", "Score", "Strand", "Sequence", "Mismatch_Number", "Mismatch_Positions","Variants")
    }
    else
    {
        colnames(offtargetData) <- c("#Chr", "Start", "End", "Name", "Score", "Strand", "Sequence", "Mismatch_Number", "Mismatch_Positions")
    }

    offtargetData[,"Score"] <- rfPrediction
    write.table(offtargetData, file=offtargetPath, quote=FALSE, sep="\t", row.names=FALSE)
    return(rfPrediction)
}

rfPredictions <- rfClassification(args[1], args[2], args[3])

