
library(caret)
resultsDir = "X:/myGit/myKaggle/Titanic/pyScripts/"

tempReal <- read.csv(paste (resultsDir, "myfirstforest.csv", sep=""))
str(tempReal)

as.REAL <- as.factor(tempReal$Survived)
str(as.REAL)


pred <- read.csv(paste (resultsDir, "genderclassmodel.csv", sep=""))
as.pred <- as.factor(pred$Survived)


pred <- read.csv(paste (resultsDir, "genderBase_submission.csv", sep=""))
as.pred <- as.factor(pred$Survived)


confusionMatrix(as.pred, as.REAL)


resultsDir = "X:/myGit/myKaggle/Titanic/rScripts/"
pred <- read.csv(paste (resultsDir, "submission.csv", sep=""))
as.pred <- as.factor(pred$Survived)


resultsDir = "X:/myGit/myKaggle/Titanic/rScripts/"
pred <- read.csv(paste (resultsDir, "1_random_forest_r_submission.csv", sep=""))
as.pred <- as.factor(pred$Survived)

