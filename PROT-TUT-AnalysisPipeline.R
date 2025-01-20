# https://rformassspectrometry.github.io/book/sec-quant.html#analysis-pipeline

library(tidyverse)
library(ggplot2)
library(QFeatures)
library(limma)

## missing values - technical or biological nature
cptac_se<-zeroIsNA(cptac_se)
nNA(cptac_se)
barplot(nNA(cptac_se)$nNAcols$pNA)
table(nNA(cptac_se)$nNArows$pNA)
cptac_se<-filterNA(cptac_se, pNA=4/6)

## missing values MAR and MNAR (random/not random) need different imputation methods

