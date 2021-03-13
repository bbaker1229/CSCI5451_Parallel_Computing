#!/usr/bin/env Rscript 

data <- read.csv("/export/scratch/users/csci5451/pollution_Vsmall.csv")
data <- data[,-c(1)]
mod <- kmeans(data, 5)
print(summary(as.factor(mod$cluster)))
