library(limma)

head(expressions)
meta<-pData(eset)

design<-model.matrix(~ 0 + meta$`disease state:ch1`)
levels(as.factor(make.names(meta$`disease state:ch1`)))
colnames(design)<- meta$`disease state:ch1` %>% make.names %>% as.factor %>% levels

fit<-lmFit(expressions, design)
contr_mat<-makeContrasts(healthy_vs_disease=multiple.sclerosis-healthy, levels=design)
fit2<-contrasts.fit(fit, contr_mat)
fit2<-eBayes(fit2)

results<-topTable(fit2, adjust.method="BH", number=Inf)
head(results)

filtered_results<-results[results$adj.P.Val<=0.05 & abs(results$logFC >= 1), ]
nrow(filtered_results)

ggplot(results, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=-log10(0.05),color="blue", linetype="dashed") +
  geom_vline(xintercept=c(-1.0,1.0), color="green", linetype="dashed") +
  labs(title="Expressions from DEA", x="Log2FC", y="-log10(adj.p.val)") +
  theme_minimal()
limma::plotMA(fit2, main="MA Plot", ylim=c(-5,5))
abline(h = c(-1, 1), col = "blue", lty = 2)

View(results)
