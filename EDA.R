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

expressions_sample_rows<-t(expressions)
expressions_umap<-umap(expressions_sample_rows)
expressions_umap_df<-as.data.frame(expressions_umap)
expressions_umap_df$Group<-pData(eset)$`disease state:ch1`
ggplot(expressions_umap_df, aes(x=V1, y=V2, colour=Group))+geom_point()

dist_expr<-dist(expressions_umap_df[,1:2], method='euclidean')
hclust_expr<-hclust(dist_expr, method="complete")
dend<-as.dendrogram(hclust_expr)
dend<-color_branches(dend, k=3)

expressions_sample_rows_df<-t(expressions)
group<-pData(eset)$`disease state:ch1`
dist_expr_non_reduced<-dist(expressions_sample_rows_df)
hclust_expr_non_reduced<-hclust(dist_expr_non_reduced, method="complete")                                                          
dend_non_reduced<-as.dendrogram(hclust_expr_non_reduced)      
dend_non_reduced<-color_branches(dend_non_reduced, k=4)

par(mfrow=c(2,1))
dend_non_reduced %>% set("labels", group) %>% set("labels_colors", as.numeric(as.factor(group)), order_value=TRUE) %>% set("labels_cex", 0.5) %>% plot()
dend %>% set("labels", expressions_umap_df$Group) %>% set("labels_colors", as.numeric(as.factor(expressions_umap_df$Group)), order_value=TRUE) %>% set("labels_cex", 0.5) %>% plot()
