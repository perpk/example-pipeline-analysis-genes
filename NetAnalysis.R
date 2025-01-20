library(WGCNA)

head(expressions)
head(deg_df)
nrow(deg_df)
expressions_deg<-expressions[rownames(expressions) %in% rownames(deg_df), ]
nrow(expressions_deg)

powers = c(1:30)
sft<-pickSoftThreshold(t(expressions_deg), powerVector=powers, verbose=5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], type="n",
     xlab="Soft Threshold (power)", ylab="Scale Free Technology Model Fit")

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.75, col="blue")
abline(h=0.8, col="green")

softPower=23
expressions_deg_t<-t(expressions_deg)
head(expressions_deg_t)
adjacency<-adjacency(expressions_deg_t, power=softPower)
TOM<-TOMsimilarity(adjacency)
dissTOM<-1-TOM
geneTree<-hclust(as.dist(dissTOM), method="average")
head(geneTree)
plot(geneTree, xlab="", sub="", main="Gene Dissimilarity Based on TOM-based Dissimilarity", labels=FALSE, hang=0.04)

dynamicMods<-cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit = 2, pamRespectsDendro=FALSE, minClusterSize=30)
table(dynamicMods)

dynamicColors<-labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuid=TRUE, guideHang=0.05)

MEList<-moduleEigengenes(expressions_deg_t, colors=dynamicColors)
MEs<-MEList$eigengenes
head(MEList)
head(MEs)

moduleMembership<-cor(expressions_deg_t, MEs, use="pairwise.complete.obs")
moduleColors<-dynamicColors

library(stringr)
library(hash)
library(writexl)
library(hgu133plus2.db)
library(AnnotationDbi)

source("./MapAffyProbeToSymbol.R")

hubGenes <- c()
h <- hash()

dim(moduleMembership)
affy_probe_to_symbol_map<-MapAffyProbeToSymbol(rownames(moduleMembership))
rownames(moduleMembership)<-affy_probe_to_symbol_map$geneSymbol

hub_genes<-data.frame()
for (group in colnames(moduleMembership)) {
  modGenes = (moduleColors == group %>% str_replace("ME",""))
  hubGene = names(which.max(moduleMembership[modGenes, group]))
  genes_in_group <- names(moduleMembership[modGenes, group])
  cat("Length for group", group, "=", length(genes_in_group), "__ ")
  h[[group]] <- names(moduleMembership[modGenes, group])
  cat("Hub gene for module", group, ":", hubGene, "\n")
  hubGenes <- c(hubGenes, hubGene)
}

affy_probe_to_symbol_map<-MapAffyProbeToSymbol(rownames(expressions_deg))
rownames(expressions_deg)<-affy_probe_to_symbol_map$geneSymbol
module="turquoise"
TOMExpr<-TOMsimilarityFromExpr(expressions_deg_t)
moduleGenes = (moduleColors == module)
expressions_deg_mod<-expressions_deg[moduleGenes == TRUE, ]
expressions_deg_t<-t(expressions_deg_mod)
cyt<-exportNetworkToCytoscape(TOM[moduleGenes, moduleGenes], 
                              edgeFile="cytoInput-edges-turquoise.txt",
                              nodeFile="cytoInput-nodes-turquoise.txt",
                              weighted=TRUE,
                              threshold=0.02,
                              nodeNames=colnames(expressions_deg_t),
                              nodeAttr = moduleColors[moduleGenes])



