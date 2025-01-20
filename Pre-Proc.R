library(GEOquery)
library(affy)
library(ggplot2)
library(car)
library(reshape2)        
library(pheatmap)
library(ggrepel)
library(tidyverse)

source("./DeleteFiles.R")

temp_d<-tempdir()
DeleteFiles(tempdir(), c("\\.CEL$", "\\.CEL.gz$"))

geods<-"GSE21942"
tar<-paste(geods, "_RAW.tar", sep="")

untar(tar, exdir=temp_d)
gz_f<-list.files(temp_d, pattern="\\CEL.gz$", full.names=TRUE)
cel_fs<-sapply(gz_f, function(f) {
  R.utils::gunzip(f, remove=FALSE)
})
celfiles<-list.files(temp_d, pattern="\\.CEL$")
affy_data<-ReadAffy(filenames=celfiles, celfile.path=temp_d)
norm_d<-rma(affy_data)
expressions<-exprs(norm_d)

# raw_vals<-exprs(affy_data)

boxplot(expressions)
expr_melt<-melt(expressions)
ggplot(expr_melt,aes(x=value, color=Var2)) +
  geom_density() +
  theme_minimal()
hist(expressions)

# raw_melt<-melt(raw_vals)
# ggplot(raw_melt,aes(x=value, color=Var2)) +
#  geom_density() +
#  theme_minimal()

gse<-getGEO(geods)
eset<-gse[[1]]
pca<-prcomp(t(expressions), scale.=TRUE)
pca_df<-as.data.frame(pca$x)
pca_df$Group<-pData(eset)$`disease state:ch1`
ggplot(pca_df, aes(x=PC1, y=PC2, color=Group, label=rownames(pca_df))) +
  geom_point(size=3) +
  geom_text_repel(size=3) +
  theme_minimal() 
  
summary(pca)
pca$sdev

var_explain<-(pca$sdev)^2 / sum((pca$sdev)^2)
cumsum_var<-cumsum(var_explain)
plot(var_explain, type='b', xlab="PC", ylab="Variance Explained", main="Scree Plot")
plot(cumsum_var, type='b', xlab="PC", ylab="Cumulative Variance Explained", main="Scree Plot - Cumulative Variance")

annotation<-data.frame(DiseaseStatus=pData(eset)$`disease state:ch1`)
rownames(annotation)<-colnames(expressions)

cor_mat<-cor(expressions)
pheatmap(cor_mat, annotation_col = annotation, 
         clustering_distance_cols = "euclidean", clustering_distance_rows = "euclidean", clustering_method="complete",
         color=colorRampPalette(c("red", "green", "blue"))(50), main="Sample correlation heatmap")



