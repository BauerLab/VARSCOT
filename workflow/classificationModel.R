########################################################################################################################
########################################################################################################################
# This notebook contains description and code of the feature matrix extraction, model building and feature selection
# based on the GUIDE-seq data by Tsai, S.Q. et al. (http://www.nature.com/nbt/journal/v33/n2/full/nbt.3117.html).

########################################################################################################################
# Data:
########################################################################################################################

library(randomForest)
library(caret)
library(ROCR)
library(ggplot2)
load(file="data-objects/datasetsSampling.RData")              # 10 datasets with downsampled inactive off-targets
source("evalFunctions.R")

########################################################################################################################
# Build feature matrix for each dataset:
########################################################################################################################

# featureMatrix <- lapply(datasetsSampling, function(x) buildFeatureMatrix(x, 23, numeric=TRUE))
# save(featureMatrix, file="data-objects/featureMatrix.RData")
load(file="data-objects/featureMatrix.RData")


########################################################################################################################
# Classification with random forest
########################################################################################################################

# Perform feature selection with all datasets combined
# Feature selection sets own seed

# rfFeatures <- rfFeatureSelectionMultiple(featureMatrix)
# save(rfFeatures, file="data-objects/rfFeatures.RData")
load(file="data-objects/rfFeatures.RData")

maxIndex <- which(rfFeatures==max(rfFeatures))
# rfSelectedFeatures <- names(rfFeatures)[maxIndex:length(rfFeatures)]
# save(rfSelectedFeatures, file="data-objects/rfSelectedFeatures.RData")
load(file="data-objects/rfSelectedFeatures.RData")

# Plot accuracies (for thesis)
png("../doc/thesis/images/feature-selection.png", width=1000, height=700)
plot(rev(rfFeatures), ylab="Accuracy based on OOB error", xlab="Feature importance rank", cex.axis=1.4, cex.lab=1.4)
abline(v=(length(rfFeatures) - maxIndex + 1), lty=2, col="Brown4", lwd=2)
dev.off()

# Final models
# set.seed(42)

# rfAllClassifiers <- list()
# for (i in 1:10)
# {
#     trainingData <- featureMatrix[[i]][,c("offtargetActivity", rfSelectedFeatures)]
#     trainingData[,"offtargetActivity"] <- as.factor(trainingData[,"offtargetActivity"])
#     rfClassifier <- randomForest(offtargetActivity~., data=trainingData, ntree=1000)
#     rfAllClassifiers <- c(rfAllClassifiers, list(rfClassifier))
# }
# save(rfAllClassifiers, file="data-objects/rfAllClassifiers.RData")
load(file="data-objects/rfAllClassifiers.RData")

################################################################################
# Figure 2a
################################################################################


# Collect average feature importance
# set.seed(42)
# variantImportanceAll <- lapply(featureMatrix, function(x)
# {
#   variantImportanceAll <- apply(do.call(cbind, lapply(1:10, function(y) randomForest(offtargetActivity~., data=x)$importance)), 1, mean)
# })
# save(variantImportanceAll, file="data-objects/variantImportanceAll.RData")
load(file="data-objects/variantImportanceAll.RData")

varImp.df = data.frame(variantImportanceAll)
colnames(varImp.df) = c(1:10)
varImp.df$Feature = row.names(varImp.df)
varImp.df$Average = rowMeans(varImp.df[,1:10], na.rm=TRUE)
for (i in c(1:nrow(varImp.df))) {
  varImp.df[i,'STDEV']   = sd((varImp.df[i,1:10]))
}

varImp.df = varImp.df[order(-varImp.df$Average),]

# select only the top 30 features
varImp.top = varImp.df[1:30,]
varImp.top$Feature = factor(varImp.top$Feature, levels=varImp.top$Feature[order(varImp.top$Average)],)

fig2a = ggplot(varImp.top, aes(x=Feature, y=Average, ordered=TRUE)) +
  geom_point() +
  ylim(limits=c(0,25)) +
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV)) +
  coord_flip()

png("images/feature-importance.png", width=600, height=500)
ggplot(varImp.top, aes(x=Feature, y=Average, ordered=TRUE)) +
  geom_point() +
  ylim(limits=c(0,25)) +
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV)) +
  coord_flip() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  labs(y = "Average mean decrease in accuracy")
dev.off()

png("../doc/thesis/images/feature-importance.png", width=1000, height=800)
ggplot(varImp.top, aes(x=Feature, y=Average, ordered=TRUE)) +
  geom_point() +
  ylim(limits=c(0,25)) +
  geom_errorbar(aes(ymin=Average-STDEV, ymax=Average+STDEV)) +
  coord_flip() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  labs(y = "Average mean decrease in accuracy in % (range for all datasets)")
dev.off()

################################################################################
# Figure 2b
################################################################################


featureMatrixTest <- featureMatrix[[1]] # just to set up variable
folds <- createFolds(1:nrow(featureMatrixTest), k = 10, list = TRUE, returnTrain = FALSE) # sets up splits for cross validation

impROCR <- list()
allROCR <- list()

for (j in 1:length(featureMatrix)) {
  print(j)
  featureMatrixTest <- featureMatrix[[j]]

  predAll <- c()
  predImp <- c()
  labels  <- c()

  a = 1

  for (i in folds) {
    print(a)
    a = a + 1
    trainingSet <- featureMatrixTest[-i,]
    testSet <- featureMatrixTest[i,]

    class = as.numeric(as.character(testSet$offtargetActivity))

    rf.imp <- randomForest(offtargetActivity~., data=trainingSet[,c("offtargetActivity", rfSelectedFeatures)], ntree=1000)
    rf.all <- randomForest(offtargetActivity~., data=trainingSet, ntree=1000)

    pred.i = predict(rf.imp, testSet, type='prob')[,2]
    pred.a = predict(rf.all, testSet, type='prob')[,2]

    predAll <- c(predAll, pred.a)
    predImp <- c(predImp, pred.i)
    labels  <- c(labels,  class)
  }

  impROCR$predictions[[j]] = predImp
  impROCR$labels[[j]]      = labels

  allROCR$predictions[[j]] = predAll
  allROCR$labels[[j]]      = labels
}

predImp <- prediction(impROCR$predictions, impROCR$labels)
predAll <- prediction(allROCR$predictions, allROCR$labels)

perfImp <- performance(predImp, 'tpr','fpr')
perfAll <- performance(predAll, 'tpr','fpr')

aucImp <- performance(predImp, measure='auc', x.measure='cutoff')
aucAll <- performance(predAll, measure='auc', x.measure='cutoff')

t.test(unlist(aucImp@y.values),unlist(aucAll@y.values),paired=TRUE)

png("images/cross-validation.png", width=700, height=500)
plot(perfAll, col='Brown3',  lty=3, lwd=1, cex.axis=1.4, cex.lab=1.4)
plot(perfAll, avg='vertical', col='Brown3', lwd=3, cex.axis=1.4, add=TRUE)
plot(perfImp, col='DeepSkyBlue3', lty=3, lwd=1, cex.axis=1.4, cex.lab=1.4, add=TRUE)
plot(perfImp, avg='vertical', col='DeepSkyBlue3', lwd=3, cex.axis=1.4, add=TRUE)
legend("bottomright", c("Selected features - single dataset", "Selected features - average datasets", "All features - single dataset", "All features - average datasets"),
col=c("DeepSkyBlue3", "DeepSkyBlue3", "Brown3", "Brown3"), lty=c(3, 1, 3, 1), cex=1.4, lwd=2)
dev.off()

png("../doc/thesis/images/cross-validation.png", width=1000, height=700)
plot(perfImp, col='DeepSkyBlue3', lty=3, lwd=1, cex.axis=1.4, cex.lab=1.4)
plot(perfImp, avg='vertical', col='DeepSkyBlue3', lwd=3, cex.axis=1.4, add=TRUE)
plot(perfAll, col='Brown3',  lty=3, lwd=1, cex.axis=1.4, add=TRUE)
plot(perfAll, avg='vertical', col='Brown3', lwd=3, cex.axis=1.4, add=TRUE)
legend("bottomright", c("Selected features - single dataset", "Selected features - average datasets", "All features - single dataset", "All features - average datasets"),
col=c("DeepSkyBlue3", "DeepSkyBlue3", "Brown3", "Brown3"), lty=c(3, 1, 3, 1), cex=1.4, lwd=2)
dev.off()

################################################################################
# Choose best classifier
################################################################################

indexBestModel <- which(unlist(aucImp@y.values)==max(unlist(aucImp@y.values)))

rfClassifier <- rfAllClassifiers[[indexBestModel]]
# save(rfClassifier, file="data-objects/rfClassifier.RData")


################################################################################
# Save best feature matrix and selected features for python script input
################################################################################

# write.table(featureMatrix[[indexBestModel]], file = "../pipeline/python-model-implementation/featureMatrixTraining.txt", quote = FALSE, sep = "\t")
# write.table(rev(rfSelectedFeatures), file = "../pipeline/python-model-implementation/selectedFeaturesRF.txt", quote = FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
