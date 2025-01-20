BiocManager::install("tidyverse")
BiocManager::install("factoextra")
BiocManager::install("msdata")
BiocManager::install("mzR")
BiocManager::install("rhdf5")#
BiocManager::install("rpx")#
BiocManager::install("MsCoreUtils")#
BiocManager::install("QFeatures")#
BiocManager::install("Spectra")#
BiocManager::install("ProtGenerics")#
BiocManager::install("PSMatch")#
BiocManager::install("pheatmap")
BiocManager::install("limma")
BiocManager::install("MSnID")
BiocManager::install("RforMassSpectrometry/SpectraVis")

library(rpx)
px <- PXDataset("PXD000001")
pxfiles(px)
pxtax(px)
pxurl(px)
pxref(px)

# The mzML format is an open, XML-based format for mass spectrometer output files, 
# developed with the full participation of vendors and researchers in order to create
# a single open format that would be supported by all software.

fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"
mzf <- pxget(px, fn)
px <- PXDataset("PXD022816")
pxget(px, grep("mzID", pxfiles(px))[1:3])
pxget(px, grep("mzML", pxfiles(px))[1:3])

BiocManager::install(c("bookdown", "BiocStyle", "lgatto/msmbstyle"))
source("https://raw.githubusercontent.com/rformassspectrometry/book/main/install_docs_deps.R")

###

library(msdata)
proteomics()
ident()
quant()

library(Spectra)
### Here, example data created on-the-fly to demonstrate a few of the Spectra capabilities.
spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
spd
sp0<-Spectra(spd)
sp0
spectraVariables(sp0)
spectraData(sp0)
peaksData(sp0)
###
sp<-Spectra(mzf)
sp
length(sp)
spectraVariables(sp)
spectraData(sp)
pd<-peaksData(sp)
View(pd)

help(DataFrame)

with(spectraData(filterMsLevel(sp, 1)), 
     plot(rtime, totIonCurrent, type="l"))
abline(v=rtime(sp)[2807], col="red")

with(spectraData(filterPrecursorScan(sp)),
     plot(MS1, MS2))

ms_2<-filterPrecursorScan(sp, 2807)
ms_2

plotSpectra(sp[2807], xlim=c(400,1000))
abline(v=precursorMz(ms_2)[-1], col="grey")
precursorMz(ms_2)[-1]
abline(v=precursorMz(ms_2)[2], col="red")
precursorMz(ms_2)[2]

plotSpectra(sp[2807], xlim=c(521.2, 522.5), type="l")
abline(v=precursorMz(ms_2)[2], col="red")

###
factor(sp$msLevel)
length(sp[msLevel=1])
length(sp[msLevel=2])

###
library(QFeatures)
data(feat1)

colData(feat1)
feat1[[1]]
feat1[["psms"]]
assay(feat1[[1]])
rowData(feat1[[1]])

feat1<-aggregateFeatures(feat1, i="psms", fcol="Sequence", name="peptides", fun=colMeans)
feat1
colMeans(assay(feat1[[1]])[1:3, ])
assay(feat1[[2]])["SYGFNAAR", ]

feat1<-aggregateFeatures(feat1, i="peptides", fcol="Protein", name="proteins", fun=colMedians)
feat1

assay(feat1[["proteins"]])
feat1
feat1["ProtA", ,]
filterFeatures(feat1, ~ pval < 0.05)

data(hlpsms)
hl<-readQFeatures(hlpsms, ecol=1:10, name="psms")
hl[[1]]
head(assay(hl[["psms"]]))
head(rowData(hl[["psms"]]))
###

basename(f<-msdata::quant(pattern="cptac", full.names=TRUE))
names(read.delim(f))

(i<-grep("Intensity\\.", names(read.delim(f))))
cptac_se<-readSummarizedExperiment(f, ecol=i, fnames="Sequence", sep="\t")
colData(cptac_se)
colnames(cptac_se)<-sub("I.+\\.", "", colnames(cptac_se))
cptac_se$condition<-sub("_[7-9]", "", colnames(cptac_se))
cptac_se$id<-sub("^.+_", "", colnames(cptac_se))
colData(cptac_se)
cptac_se

keepVar<-c("Sequence", "Proteins", "Leading.razor.protein", "PEP", "Score",
           "Reverse", "Potential.contaminant")
rowData(cptac_se)<-rowData(cptac_se)[, keepVar]
rowData(cptac_se)

