
setwd("/Users/li11/myGit/workingWithYicheng/phase-I-data/reconData")
dt <- read.table (f[1], header = T, sep = "\t")

set.seed (123456)
select.row <- c(sample (which(dt$Class == 1), 40),  sample (which(dt$Class == 0), 40))

train <- dt[select.row, ]
write.table (train, "training-data.txt", sep = "\t")

test <- dt [-select.row,]
write.table (test, "testing-data.txt", sep = "\t")

test$Class = 1
write.table (test, "data-2-test.txt", sep = "\t")

