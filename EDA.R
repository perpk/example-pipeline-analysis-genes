# hierarchical clusering

View(expressions)

dist<-as.dist(1-cor(expressions))
hcg<-hclust(dist, method="average")
plot(hcg)

# UMAP or t-SNE
install.packages("uwot")
install.packages("dendextend")
library(umap)
library(uwot)
library(ggplot2)
library(dendextend)
library(factoextra)
library(cluster)
library(tidyverse)

expressions_sample_rows<-t(expressions)
expressions_umap<-umap(expressions_sample_rows)
expressions_umap_df<-as.data.frame(expressions_umap)
expressions_umap_df$Group<-pData(eset)$`disease state:ch1`
ggplot(expressions_umap_df, aes(x=V1, y=V2, colour=Group))+geom_point()

dist_expr<-dist(expressions_umap_df[,1:2], method='euclidean')
hclust_expr<-hclust(dist_expr, method="complete")
dend<-as.dendrogram(hclust_expr)
dend<-color_branches(dend, k=2)

expressions_sample_rows_df<-data.frame(t(expressions))
#group<-pData(eset)$`disease state:ch1`
#expressions_sample_rows_df<-cbind(expressions_sample_rows_df, pData(eset)$`disease state:ch1`)
expressions_sample_rows_df$Group<-pData(eset)$`disease state:ch1`
dist_expr_non_reduced<-dist(expressions_sample_rows_df[, !names(expressions_sample_rows_df) %in% "Group"], method="euclidean")
hclust_expr_non_reduced<-hclust(dist_expr_non_reduced, method="complete")                                                          
dend_non_reduced<-as.dendrogram(hclust_expr_non_reduced)      
dend_non_reduced<-color_branches(dend_non_reduced, k=2)

#par(mar=c(1,1,1,1))
#par(mfrow=c(2,1))
dend_non_reduced %>% set("labels", expressions_sample_rows_df$Group) %>% 
  set("labels_colors", as.numeric(as.factor(expressions_sample_rows_df$Group)), order_value=TRUE) %>% 
  set("labels_cex", 0.5) %>% 
  plot()
dend %>% set("labels", expressions_umap_df$Group) %>% 
  set("labels_colors", as.numeric(as.factor(expressions_umap_df$Group)), order_value=TRUE) %>% 
  set("labels_cex", 0.5) %>% 
  plot()

fviz_silhouette(silhouette(cutree(hclust_expr, 2), dist_expr))
fviz_silhouette(silhouette(cutree(hclust_expr, 3), dist_expr))
fviz_silhouette(silhouette(cutree(hclust_expr, 4), dist_expr))
fviz_silhouette(silhouette(cutree(hclust_expr, 5), dist_expr))

fviz_silhouette(silhouette(cutree(hclust_expr_non_reduced, 2), dist_expr_non_reduced))
fviz_silhouette(silhouette(cutree(hclust_expr_non_reduced, 3), dist_expr_non_reduced))
fviz_silhouette(silhouette(cutree(hclust_expr_non_reduced, 4), dist_expr_non_reduced))
fviz_silhouette(silhouette(cutree(hclust_expr_non_reduced, 5), dist_expr_non_reduced))
