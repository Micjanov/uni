###################################
# Practicum 2: Microarrays Part 2 #
###################################



############
# Infinium #
############


## Set up your working directory
#####################

setwd("C:/Users/Louis Coussement/Documents/AAP/Education/Teaching/AHTA/B_AppliedHighthroughputAnalysis/Datasets/Prac_02/Infinium")

## Instal packages
#####################

## Install packages
BiocManager::install('lumi',update=F)
BiocManager::install('wateRmelon',update=F)
BiocManager::install('minfi',update=F)
BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b2.hg19',update=F)
BiocManager::install('IlluminaHumanMethylationEPICmanifest',update=F)
BiocManager::install('ChAMPdata',update=F)


## Load packages
library('lumi')
library('wateRmelon')
library('ChAMPdata')


## Load and explore dataset
######################

## Load the Infinium data
infdata_all <- readEPIC(getwd())

## Load the annotation data in:
annot <- read.table("Annotation_Infinium.txt",header=T,sep="\t")

## Have a look at the data and annotation
print(infdata_all)
print(dim(infdata_all))
print(annot)
print(sum(is.na(exprs(infdata_all))))
print(head(betas(infdata_all)))
print(head(exprs(infdata_all)))


## Switch to output directory
#####################

setwd("C:/Users/Louis Coussement/Documents/AAP/Education/Teaching/AHTA/B_AppliedHighthroughputAnalysis/Output/Prac_02")


## Preprocessing 
#######################

## Remove NA values
infdata_all <- infdata_all[rowSums(is.na(exprs(infdata_all)))==0,]

## Take subset_al
infdata <- infdata_all[,c(grep("CAF", annot[,2]), grep("NAF", annot[,2]))]
annot <- annot[c(grep("CAF",annot[,2]),grep("NAF",annot[,2])),]

## Change sampleNames to something more comprehensible
sampleNames(infdata) <- paste(annot[,1],annot[,2],sep="_")

## Remove probes for which calling p-value is insufficient
infdata.pf<-pfilter(infdata)

## Comparison of average methylation between control and tumor samples
# Low dimensional data
jpeg("OverallMethylation.jpg")
boxplot(betas(infdata),las=2)
dev.off()
# High dimensional data 
meth_mean_CAF <- rep(0,3)
meth_mean_NAF <- rep(0,3)
for (i in 1:ncol(infdata)){
  if (i < 4){
    meth_mean_CAF[i] <- mean(betas(infdata)[,i])
  } else {
    meth_mean_NAF[i-3] <- mean(betas(infdata)[,i])
  }
}

t_test_res <- t.test(meth_mean_CAF,meth_mean_NAF)
t_test_res
dat_boxplot <- data.frame(betas = c(meth_mean_CAF,meth_mean_NAF),
                          group = c("CAF","CAF","CAF","NAF","NAF","NAF"))
jpeg("OverallMethylation_averages.jpg")
par(mfrow=c(1,2))
boxplot(betas~group,dat_boxplot,las=2)
boxplot(betas~group,dat_boxplot,las=2,ylim=c(0,1))
dev.off()


## Normalization & QC
#######################

## Perform normalization including dye color adjustment
infdata.dasen.pf <- dasen(infdata.pf)

## Make methylumi objects to check density and color bias adjustment
infdataM <- as(infdata.pf, 'MethyLumiM')
infdataN <- as(infdata.dasen.pf, 'MethyLumiM')

## Make QC plot
jpeg("QC.jpg")
par(mfrow=c(2,2))
plotColorBias1D(infdataM,channel="both",main="before")
plotColorBias1D(infdataN,channel="both",main="after")
density(infdataM,xlab="M-value",main="before")
density(infdataN,xlab="M-value",main="after")
dev.off()


## Differential methylation analysis: limma
############################

## Build design and contrasts
des <- factor(as.character(annot[,2]))
design <- model.matrix(~0+des)
colnames(design) <- c("Tumor","Normal")
cont.matrix <- makeContrasts(NvsS=Tumor-Normal,levels=design)

## Limma
fit <- lmFit(infdataN,design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)

## Volcano plot
jpeg("Volcanoplot.jpg")
volcanoplot(fit2)
dev.off()

## Volcano plot
jpeg("MAplot.jpg")
limma::plotMA(fit2)
dev.off()

# DE results
LIMMAout <- topTable(fit2,adjust="BH",number=nrow(exprs(infdataM)))
head(LIMMAout)

## Check M-values for top results
exprs(infdataN)[rownames(infdataN)%in%rownames(head(LIMMAout)),]
# note that the probesets are not necessarly in the same order as for the limma output


## Functional annotation of limma results
############################

## Load annotation and sort alphabetically on probe name
data("probe.features.epic")
annotation_MA <- probe.features
print(head(annotation_MA))
annotation_MA <- annotation_MA[sort(rownames(annotation_MA),index.return=T)$ix,]

## Check if all probes are present in both sets
dim(LIMMAout)
dim(annotation_MA)
sum(LIMMAout$Probe_ID%in%rownames(annotation_MA))
sum(rownames(annotation_MA)%in%LIMMAout$Probe_ID) 
# Also check the reverse so no duplicate rows are present in annotation

## Since more probes are present in the annotation file, remove unnecessary probes
annotation_MA <- annotation_MA[rownames(annotation_MA)%in%LIMMAout$Probe_ID,]

## Sort LIMMA output alphabetically on probe name
LIMMAout_sorted <- LIMMAout[sort(LIMMAout$Probe_ID,index.return=T)$ix,]

## Add gene names to LIMMA output
LIMMAout_sorted$Gene <- annotation_MA$gene
LIMMAout_sorted$Feature <- annotation_MA$feature
LIMMAout_sorted$Chrom <- annotation_MA$CHR
LIMMAout_sorted$Pos <- annotation_MA$MAPINFO
LIMMAout_sorted$Chrom <- as.character(LIMMAout_sorted$Chrom)
LIMMAout_sorted$Gene <- as.character(LIMMAout_sorted$Gene)
LIMMAout_sorted$Feature <- as.character(LIMMAout_sorted$Feature)
# The data type for these columns is altered to prevent issues further downstream


## Quantification of absolute methylation differences
############################

## Check if dimension of objects to combine are the same
dim(LIMMAout_sorted)
dim(betas(infdata))

## Add gene names to LIMMA output
LIMMAout_sorted$Tumor_meth <- rowMeans(betas(infdata)[rownames(infdata)%in%LIMMAout_sorted$Probe_ID,annot$Origin=="CAF"])
LIMMAout_sorted$Control_meth <- rowMeans(betas(infdata)[rownames(infdata)%in%LIMMAout_sorted$Probe_ID,annot$Origin=="NAF"])
LIMMAout_sorted$Abs_diff_meth <- abs(rowMeans(betas(infdata)[rownames(infdata)%in%LIMMAout_sorted$Probe_ID,annot$Origin=="CAF"]) -
                                       rowMeans(betas(infdata)[rownames(infdata)%in%LIMMAout_sorted$Probe_ID,annot$Origin=="NAF"]))


## Resort results
############################

LIMMAout_annot <- LIMMAout_sorted[sort(LIMMAout_sorted$P.Value,index.return=T)$ix,c(1,12,13,10,11,4,7,8,5,14,15,16)] 
# Sort on p-values to prevent errors in sorting due to equal FDR values


## Interpretation results
############################

## Select CpGs in genic regions
sum(LIMMAout_annot$adj.P.Val<0.05)
sum(LIMMAout_annot$adj.P.Val[LIMMAout_annot$Gene!=""]<0.05)

LIMMAout_annot_gene <- LIMMAout_annot[LIMMAout_annot$Gene!="",]

## Check genic results 
head(LIMMAout_annot_gene[c(4,5,6,8,10,11,12)])
# these columns are selected because they contain the most interesting fields, however you can of course return all as well
topgenes_genic <- unique(LIMMAout_annot_gene$Gene[1:10])
for (i in 1:length(topgenes_genic)){
  LIMMAout_subset <- LIMMAout_annot_gene[(LIMMAout_annot_gene$Gene==topgenes_genic[i]) & 
                                           (LIMMAout_annot_gene$adj.P.Val<0.05) & (abs(LIMMAout_annot_gene$logFC)>2),]
  print(LIMMAout_subset[sort(LIMMAout_subset$Pos,index.return=T)$ix,c(4,5,6,8,10,11,12)])
}

## Select CpGs in promoter regions
LIMMAout_annot_prom <- LIMMAout_annot_gene[grepl("TSS",LIMMAout_annot_gene$Feature) | (LIMMAout_annot_gene$Feature=="1stExon"),]
head(LIMMAout_annot_prom)

## Look for multiple CpG in promoter regions undergoing similar methylation differences
topgenes_prom <- unique(LIMMAout_annot_prom$Gene[1:10])
for (i in 1:length(topgenes_prom)){
  LIMMAout_subset <- LIMMAout_annot_prom[(LIMMAout_annot_prom$Gene==topgenes_prom[i]) & 
                                           (LIMMAout_annot_prom$adj.P.Val<0.10),]
  if(nrow(LIMMAout_subset)>1){
    print(LIMMAout_subset[sort(LIMMAout_subset$Pos,index.return=T)$ix,c(4,5,6,8,10,11,12)])
  }
}




#########################
# SNP-array  (optional) #
#########################


## Set up your working directory
#####################

setwd("C:/Users/Louis Coussement/Documents/AAP/Education/Teaching/AHTA/B_AppliedHighthroughputAnalysis/Datasets/Prac_02/SNP_data/")


## Instal packages
#####################


## Install packages
BiocManager::install('crlmm',update=F)
BiocManager::install('hapmapsnp6',update=F)
BiocManager::install('oligoClasses',update=F)
BiocManager::install('genomewidesnp6Crlmm',update=F)

## Load packages
library('crlmm')
library('hapmapsnp6')
library('oligoClasses')
library('genomewidesnp6Crlmm')


## Load and explore dataset
######################

## Load in dataset (only male individuals)
crlmmout<-crlmm(dir(),gender=c(1,1))

## Have a look at the calls and the confidences on those calls
head(oligoClasses::calls(crlmmout)) #1 = AA, 2 = AB, 3 = BB
head(confs(crlmmout))
#oligoClasses:: is added since there is a calls function in one of the packages used in subsequent analysis

## Calculate the heterozygote fractions as indication of LOH
s1<-table(as.numeric(oligoClasses::calls(crlmmout)[,1]))
s2<-table(as.numeric(oligoClasses::calls(crlmmout)[,2]))
s1[2]/sum(s1)
s2[2]/sum(s2)

## Statistical test (chi-square test:)
chisq.test(rbind(c(s1[2],sum(s1)-s1[2]),c(s2[2],sum(s2)-s2[2])))




########
# aCGH #
########


## Set up your working directory
#####################

setwd("C:/Users/Louis Coussement/Documents/AAP/Education/Teaching/AHTA/B_AppliedHighthroughputAnalysis/Datasets/Prac_02/aCGH")


## Instal packages
#####################

## Install packages
install.packages("GEOquery", update=F)
BiocManager::install("CGHcall", update=F)
BiocManager::install("CGHregions", update=F)

## Load packages
library(GEOquery)
library(CGHcall)
library(CGHregions)


## Load and explore dataset
######################

## Load data
acgh <- getGEO('GSE54118')

## Check if acgh is list
typeof(acgh)

## If it is a list run code below to get the S4 object
acgh <- acgh[[1]]
typeof(acgh)

## Have a look at the data
head(exprs(acgh))


## Load annotation and get the right format for downstream analysis
#########################

## Load in annotation (text file can be found on Minerva)
annot <- read.table('Annotation_aCGH.txt',header=T,sep='\t')

## Get the right format
annot <- annot[,1:4]
rownames(annot) <- annot[,1]
annot[,2] <- substr(annot[,2],4,5:6)
dat_cgh <- cbind(annot[rownames(acgh),],exprs(acgh))
rownames(dat_cgh) <- NULL

## Check format
head(dat_cgh)
dim(dat_cgh)

## Select only data from chromosome 1 untill 28
dat_cgh <- dat_cgh[dat_cgh[,2]%in%c(1:28),]
dim(dat_cgh)

## make cghRaw object 
cgh <- make_cghRaw(dat_cgh)
# Remark: if NA's are produced probably no chromosome filtering


## Preprocessing, normalization and segmentation
######################

## Preprocessing
cgh_preproc <- preprocess(cgh)

## Normalization
cgh_norm <- normalize(cgh_preproc)

## Segmentation
cgh_seg <- segmentData(cgh_norm)

## Normalization
cgh_segPN <- postsegnormalize(cgh_seg)


## Calling of CNV
####################

## Calling
cgh_res <- CGHcall(cgh_segPN)
cgh_res <- ExpandCGHcall(cgh_res,cgh_segPN)


## Visual representation
#######################

jpeg("aCGH_A.jpg")
par(mfrow=c(3,1))
plot(cgh_res[,1])
plot(cgh_res[,2])
plot(cgh_res[,3])
dev.off()

regions <- CGHregions(cgh_res)
jpeg("aCGH_B.jpg")
plot(regions)
dev.off()
