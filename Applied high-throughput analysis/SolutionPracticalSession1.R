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



## Import Data using Oligo
####################

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

## Set the working directory to the directory of your choice

## Additional preprocessing
pData(RMA)
head(exprs(RMA))

targets <- read.table("targets.txt", header=T)
subject <- factor(targets$sample)
case <- factor(targets$case, levels=c("control","case"))
annot <- model.matrix(~case)

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
hist(LIMMAout$P.Value)
hist(LIMMAout$adj.P.Val)
LIMMAoutlfc <- topTable(fit2,adjust="BH",number=nrow(exprs(RMA)), lfc=1)

## Check intensity values for top results
exprs(RMA)[rownames(exprs(RMA))%in%rownames(head(LIMMAout)),]
rowMeans(exprs(RMA)[rownames(exprs(RMA))%in%rownames(head(LIMMAout)),c(1,6,7)])
rowMeans(exprs(RMA)[rownames(exprs(RMA))%in%rownames(head(LIMMAout)),2:5])

