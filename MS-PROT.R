install.packages("wrProteo")
install.packages("DEP")
library(wrProteo)

dataMQ<-readMaxQuantFile(path=".", fileName = "mq-prot-search.txt")
summary(dataMQ$quant)
dataMQ
