## Set the working directory
setwd("F:/Unief/Applied High Throughput Analysis/Project/ahta-project/array")


## Load packages
library(arrayQualityMetrics)
library(ArrayExpress)
library(limma)
library(siggenes)
library('lumi')
library('wateRmelon')
library('ChAMPdata')
library("oligo")
library(affycoretools)



## Import Data using Oligo (data already downloaded)
####################

#celpath must only contain CEL files!
celpath  = "F:/Unief/Applied High Throughput Analysis/Project/ahta-project/array/cel"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)

## Having a look at the data 
head(exprs(data))
pData(data)


## Quality Control on raw data
####################

arrayQualityMetrics(data,outdir="F:/Unief/Applied High Throughput Analysis/Project/ahta-project/array/qc/raw",force=T)
arrayQualityMetrics(data,outdir="F:/Unief/Applied High Throughput Analysis/Project/ahta-project/array/qc/rawlog", force=T,do.logtransform=T)


## Preprocessing
###################

## RMA
RMA <- oligo::rma(data,background=T)
## Annotation
RMA <- annotateEset(RMA, pd.hugene.2.1.st)


## Quality Control on preprocessed data
####################

## QC post preprocessing
arrayQualityMetrics(RMA,outdir="F:/Unief/Applied High Throughput Analysis/Project/ahta-project/array/qc/postpro",force=T)  		


## Differential expression analysis with RMA preprocessed data
####################

pData(RMA)
head(exprs(RMA))

targets <- read.table("targets.txt", header=T)
subject <- factor(targets$sample)
case <- factor(targets$case, levels=c("control","case"))

design <- model.matrix(~0+case)
colnames(design)<-c("control","case")
fit <- lmFit(RMA,design)
cont.matrix <- makeContrasts(CasevsControl=case-control,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix) 
fit2 <- eBayes(fit2)

# Volcano plot
jpeg("Volcanoplot.jpg")
volcanoplot(fit2)
dev.off()

# MA plot
jpeg("MAplot.jpg")
limma::plotMA(fit2)
dev.off()
  
# DE results
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(RMA)))
head(LIMMAout)
jpeg("pvalue.jpg")
hist(LIMMAout$P.Value)
dev.off()
jpeg("qvalue.jpg")
hist(LIMMAout$adj.P.Val)
dev.off()
LIMMAoutna <- LIMMAout[complete.cases(LIMMAout), ]
LIMMAoutp <- LIMMAout[LIMMAout$adj.P.Val < 0.05,]
jpeg("lfc.jpg")
hist(LIMMAoutp$logFC, breaks=20)
dev.off()
LIMMAoutlfc1 <- topTable(fit2,adjust="BH",number=nrow(exprs(RMA)), lfc=1)
results <- LIMMAoutlfc[complete.cases(LIMMAoutlfc), ]