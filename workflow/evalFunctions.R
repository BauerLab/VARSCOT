########################################################################################################################
# Functions:
########################################################################################################################

# Build feature matrix for a dataset either for training (datasets produced by processDataForModel.R)
# or for testing (datasets for validation). Input datasets for training need a class column in the dataset.
buildFeatureMatrix <- function(data, seqLength, numeric=FALSE, validate=FALSE)
{
    # Define DNA letters, DNA pairs and transitions/transversions.
    dnaLetters <- c("A","C","G","T")
    allDnaPairs <- sort(levels(interaction(expand.grid(dnaLetters, dnaLetters),sep="")))
    transitions <- c("AG","GA","CT","TC")
    transversions<- c("AC","CA","AT","TA","GC","CG","GT","TG")

    # Features: mismatch positions and mismatch types
    mismatchPositions <- matrix(0, nrow=nrow(data), ncol=(seqLength-2))
    mismatchTypes <- matrix(0, nrow=nrow(data), ncol=12)
    colnames(mismatchTypes) <- sort(c(transitions, transversions))

    # Features: transition and transversion number
    transitionNumber <- matrix(0, nrow=nrow(data), 1)
    transversionNumber <- matrix(0, nrow=nrow(data), 1)

    adjMismatchNumber <- rep(0, nrow(data))
    seedMismatchCount <- rep(0, nrow(data))

    for (i in 1:nrow(data))
    {
        ontarget <- strsplit(data[i,"Target_Sequence"], "")[[1]]
        offtarget <- strsplit(data[i,"Offtarget_Sequence"], "")[[1]]
        for (j in 1:(seqLength - 2))
        {
            if (ontarget[j] != offtarget[j])
            {
                mismatch <- paste(ontarget[j], offtarget[j], sep="")
                mismatchPositions[i,j] <- 1
                mismatchTypes[i,mismatch] <- 1
                if (is.element(mismatch, transitions)) transitionNumber[i,1] <- transitionNumber[i,1] + 1
                else transversionNumber[i,1] <- transversionNumber[i,1] + 1

                if (j > 8 && j < 21)
                {
                    seedMismatchCount[i] <- seedMismatchCount[i] + 1
                }
                if ((j > 1) && (ontarget[j - 1] != offtarget[j - 1]))
                {
                    adjMismatchNumber[i] <- adjMismatchNumber[i] + 1
                }
            }
        }
    }

    # Features: first PAM letter
    pamLetterStart <- matrix(0, nrow=nrow(data), ncol=4)
    colnames(pamLetterStart) <- dnaLetters

    # Features: Single base and adjacent pairs
    singleLettersOfftarget <- matrix(0, nrow=nrow(data), ncol=4*(seqLength-3))
    colnames(singleLettersOfftarget) <- do.call(c, lapply(1:(seqLength-3), function(x) paste(dnaLetters, x, sep="")))

    pairedLettersOfftarget <- matrix(0, nrow=nrow(data), ncol=16*(seqLength-4))
    colnames(pairedLettersOfftarget) <- do.call(c, lapply(1:(seqLength-4), function(x) paste(allDnaPairs, x, sep="")))

    # Features: number of adjacent pairs
    pairedLettersOfftargetNumber <- matrix(0, nrow=nrow(data), ncol=16)
    colnames(pairedLettersOfftargetNumber) <- allDnaPairs

    for (i in 1:nrow(data))
    {
        offtarget <- strsplit(data[i,"Offtarget_Sequence"], "")[[1]]
        for (j in 1:(length(ontarget)-4)) # Without PAM and last letter because of pairs (last letter is added later on)
        {
            currentLetterOfftarget <- paste(offtarget[j], j, sep="")
            singleLettersOfftarget[i,currentLetterOfftarget] <- 1

            currentPairOfftarget <- paste(offtarget[j], offtarget[j+1], sep="")
            currentPairOfftargetPos <- paste(currentPairOfftarget, j, sep="")
            pairedLettersOfftarget[i,currentPairOfftargetPos] <- 1
            pairedLettersOfftargetNumber[i,currentPairOfftarget] <- pairedLettersOfftargetNumber[i,currentPairOfftarget] + 1
        }
        currentLetterOfftarget <- paste(offtarget[j+1], j+1, sep="")
        singleLettersOfftarget[i,currentLetterOfftarget] <- 1
        pamLetterStart[i,offtarget[j+2]] <- 1
    }

    # Build complete feature matrix
    if (validate)
    {
        featureMatrix <- cbind(data$NM, mismatchPositions, mismatchTypes, transitionNumber, transversionNumber, singleLettersOfftarget,
                               pamLetterStart, pairedLettersOfftarget, pairedLettersOfftargetNumber, adjMismatchNumber, seedMismatchCount, data$Target_Activity)
        colnames(featureMatrix) <- c("totalMismatches", paste("mismatchPos", 1:(seqLength-2), sep=""),
                                    sapply(colnames(mismatchTypes), function(x) paste(strsplit(x, "")[[1]][1], "to", strsplit(x, "")[[1]][2], sep="")),
                                    "transitionNumber", "transversionNumber", colnames(singleLettersOfftarget), paste("PAM", colnames(pamLetterStart), sep=""), colnames(pairedLettersOfftarget),
                                    colnames(pairedLettersOfftargetNumber), "adjacentMismatches", "seedMismatches", "ontargetActivity")
    }
    else
    {
        featureMatrix <- cbind(data$Class, data$NM, mismatchPositions, mismatchTypes, transitionNumber, transversionNumber, singleLettersOfftarget,
                               pamLetterStart, pairedLettersOfftarget, pairedLettersOfftargetNumber, adjMismatchNumber, seedMismatchCount, data$Target_Activity)
        colnames(featureMatrix) <- c("offtargetActivity", "totalMismatches", paste("mismatchPos", 1:(seqLength-2), sep=""),
                                    sapply(colnames(mismatchTypes), function(x) paste(strsplit(x, "")[[1]][1], "to", strsplit(x, "")[[1]][2], sep="")),
                                    "transitionNumber", "transversionNumber", colnames(singleLettersOfftarget), paste("PAM", colnames(pamLetterStart), sep=""), colnames(pairedLettersOfftarget),
                                    colnames(pairedLettersOfftargetNumber), "adjacentMismatches", "seedMismatches", "ontargetActivity")
    }

    featureMatrix <- apply(featureMatrix, 2, as.numeric)
    featureMatrix <- as.data.frame(featureMatrix)

    # Numeric option for classifiers that cannot deal with factors
    if (numeric && !validate)
    {
        featureMatrix[,1] <- as.factor(featureMatrix[,1])
        return(featureMatrix)
    }
    else if (!numeric && !validate)
    {
        binaryIndex <- c(1, 4:35, 38:425)
        featureMatrix[,binaryIndex] <- lapply(featureMatrix[,binaryIndex], as.factor)
    }
    else if (!numeric && validate)
    {
        binaryIndex <- c(4:35, 38:425)
        featureMatrix[,binaryIndex] <- lapply(featureMatrix[,binaryIndex], as.factor)
    }
    return(featureMatrix)
}


## Feature selection with all datasets and RF
## Sets its own seed
rfFeatureSelectionMultiple <- function(dataList)
{
    set.seed(42)
    variantImportanceAll <- lapply(dataList, function(x)
    {
        variantImportanceAll <- apply(do.call(cbind, lapply(1:10, function(y) randomForest(offtargetActivity~., data=x)$importance)), 1, mean)
    })

    ## Mean importance over all datasets
    meanImportance <- apply(do.call(cbind, variantImportanceAll), 1, mean)
    meanImportanceNames <- names(meanImportance)[order(meanImportance)]

    rankingAllDatasets <- matrix(NA, nrow=(ncol(dataList[[1]])-1), ncol=length(dataList))
    for (j in 1:length(dataList))
    {
        print(j)
        ranking <- c()
        for (i in 1:length(meanImportance))
        {
            ranking <- c(ranking, mean(unlist(lapply(1:5, function(x) 1-sum(randomForest(offtargetActivity~., data=dataList[[j]][,c("offtargetActivity", meanImportanceNames[i:length(meanImportanceNames)])])$confusion[,"class.error"])))))
        }
        rankingAllDatasets[,j] <- ranking
    }

    meanRanking <- apply(rankingAllDatasets, 1, mean)
    names(meanRanking) <- meanImportanceNames
    return(meanRanking)
}


## Hamming distance between two character strings of equal length

stringDiff <- function(x, y)
{
    if (nchar(x)!=nchar(y))
    {
        stop("Error: Strings must have the same length.")
    }
    x_vec <- strsplit(x, "")[[1]]
    y_vec <- strsplit(y, "")[[1]]

    count = 0
    for (i in 1:length(x_vec))
    {
        if (x_vec[i]!=y_vec[i]) count = count + 1
    }
    return(count)
}

# CCTop scoring
# Sum over all mismatches, compute 1.2^pos in each step
cctopScoreFunction <- function(ontarget, offtarget)
{
    if (nchar(ontarget)!=nchar(offtarget))
    {
        stop("Error: Strings must have the same length.")
    }
    ontargetVec <- strsplit(ontarget, "")[[1]]
    offtargetVec <- strsplit(offtarget, "")[[1]]

    score <- 0
    for (i in 1:length(ontargetVec))
    {
        if (ontargetVec[i]!=offtargetVec[i])
        {
            score <- score + (1.2^i)
        }
    }
    return(score)
}
