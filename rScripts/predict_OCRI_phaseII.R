##============================================
##  File: predict_OCRI_phaseII.R
##  Author: Jianying
##  Modified from predMod_3classes.R
##============================================
library(caret)
library(pROC)
library(Metrics)


#===========================
mac.os  <- "/Users/li11/"
linux   <- "~/"
windows <- "X:/"

#root <- windows
root <- mac.os
#==========================



require(compiler)
multiClassSummary <- cmpfun(function (data, lev = NULL, model = NULL)
{
                            #Load Libraries
    require(Metrics)
    require(caret)
    #Check data
    if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
      stop("levels of observed and predicted data do not match")
      #Calculate custom one-vs-all stats for each class
     prob_stats <- lapply(levels(data[, "pred"]), function(class)
       {
          #Grab one-vs-all data for the class
           pred <- ifelse(data[, "pred"] == class, 1, 0)
           obs <- ifelse(data[, "obs"] == class, 1, 0)
           prob <- data[,class]
         #Calculate one-vs-all AUC and logLoss and return
           cap_prob <- pmin(pmax(prob, .000001), .999999)
           prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
           names(prob_stats) <- c("ROC", "logLoss")
           return(prob_stats)
      })
      prob_stats <- do.call(rbind, prob_stats)
      rownames(prob_stats) <- paste( "Class:" , levels(data[, "pred"]))
    
      #Calculate confusion matrix-based statistics
      CM <- confusionMatrix(data[, "pred"], data[, "obs"])
                            #Aggregate and average class-wise stats
                            #Todo: add weights
      class_stats <- cbind(CM$byClass, prob_stats)
      class_stats <- colMeans(class_stats)
                            #Aggregate overall stats
      overall_stats <- c(CM$overall)
    
      #Combine overall with class-wise stats and remove some stats we don't want
      stats <- c(overall_stats, class_stats)
      stats <- stats[! names(stats) %in% c("AccuracyNull","Prevalence", "Detection Prevalence")]
                            
    #Clean names and return
                  
    names(stats) <- gsub('[[:blank:]] +', '_' , names(stats))
    return(stats)
})

## Note: no visible binding for global variable 'Metrics'
## Note: no visible binding for global variable 'caret' 

## set up working directory


setwd(paste (root, "/myGit/workingWithYicheng/phase-I-data/reconData", sep=""))

#setwd(paste (root, "/myGit/workingWithYicheng/phase-II-data/reconData/", sep=""))




##	param4
data <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")
setwd(paste (root, "/myGit/workingWithYicheng/phase-II-data/modeling/", sep=""))
#sink ("log_param4.txt")                   


##	data cleaning

var0 <- unlist(lapply(data, function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN0 <- data[,-which(var0)]
# drop the first column of ID?
dataN0[,1] <- NULL
data.2.classes <- dataN0[-which (dataN0$label == "k"),]
labelTrain <- data.2.classes



##### BEGIN: data partition >>>>>
## set random seed
#set.seed(12345)


## create data partition
#inTrainingSet <- createDataPartition(data$label, p=.7, list=FALSE)

dim(labelTrain)

##### BEGIN: tune the parameters >>>>>
## control:
# resampling technique: 5-repeat 10-fold cross-validation
# performance metrics: ROC AUC curve

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     summaryFunction = multiClassSummary,
                     classProbs = TRUE)


##### END: tune the parameters <<<<<
##### BEGIN: train model - svm >>>>>

set.seed(1024)
svmFit <- train(label ~ ., data = labelTrain,
                ## training model: svm >>>
                method = "svmRadial",
                metric = "ROC",
                tuneLength = 10,
                trControl = ctrl)
##=================================================================================

setwd(paste (root, "/myGit/workingWithYicheng/phase-II-data/reconData/", sep=""))
dt.phase.II <- read.table("recon_3classes_para4.txt", header=TRUE, sep = "\t")

var1 <- unlist(lapply(dt.phase.II , function(x) 0 == var(if (is.factor(x)) as.integer(x) else x)))
dataN1 <- dt.phase.II[,-which(var1)]
# drop the first column of ID?
dataN1[,1] <- NULL

labelTest <- dataN1
dim(labelTest)

#nrow(labelTrain)
#nrow(labelTest)



##### BEGIN: tune the parameters >>>>>
## control:
# resampling technique: 5-repeat 10-fold cross-validation
# performance metrics: ROC AUC curve

ctrl <- trainControl(method = "repeatedcv",
                     repeats = 5,
                     summaryFunction = multiClassSummary,
                     classProbs = TRUE)


##### END: tune the parameters <<<<<
##### BEGIN: train model - svm >>>>>

set.seed(1024)
svmFit <- train(label ~ ., data = labelTrain,
                ## training model: svm >>>
                method = "svmRadial",
                metric = "ROC",
                tuneLength = 10,
                trControl = ctrl)
## prediction
svmPred <- predict(svmFit, labelTest)
#str(svmPred)

## predicted probabilities
svmProbs <- predict(svmFit, labelTest, type = "prob")

#apply (svmProbs, 1, sum)

#str(svmProbs)
cat ("This is the prediction with SVM")
cat("\n")
cat("\n")

lab.labels <- c (rep ("c", 42), rep ("k",27), rep("n",0))

confusionMatrix(svmPred, as.factor(lab.labels))



