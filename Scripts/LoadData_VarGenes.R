# Analysis of 10X data from Eirini-Josi collaboration
library(ggplot2)
library("Seurat")
library("reshape2")
library("pheatmap")
library("plyr")
library("dplyr")
library("colorRamps")
library("SDMTools")
library("fields")
library("magrittr")
library(beepr)
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R") - To install cellRanger kit
library(cellrangerRkit)

# 1- Load matrices

setwd("~/Desktop/AML_Eirini_Josi/")
AML_DAY1_N1 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day1_n1/")
AML_DAY1_N2 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day1_n2/")
AML_DAY7_N1 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day7_n1/")
AML_DAY7_N2 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day7_n2/")
AML_DAY14_N1 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day14_n1/")
AML_DAY14_N2 <- load_cellranger_matrix("~/Desktop/Chromium 10x data/AML-4_10day14_n2/")

# 2- Merge data into one full matrix and create annotation vector

all=cbind(exprs(AML_DAY1_N1),exprs(AML_DAY1_N2),exprs(AML_DAY7_N1),exprs(AML_DAY7_N2),exprs(AML_DAY14_N1),exprs(AML_DAY14_N2))
mapping = c(rep("A.Day1_N1", ncol(AML_DAY1_N1)), rep("B.Day1_N2", ncol(AML_DAY1_N2)), rep("C.Day7_N1", ncol(AML_DAY7_N1)), rep("D.Day7_N2", ncol(AML_DAY7_N2)), rep("E.Day14_N1", ncol(AML_DAY14_N1)), rep("F.Day14_N2", ncol(AML_DAY14_N2)))
day = c(rep("A.Day1", ncol(AML_DAY1_N1)), rep("A.Day1", ncol(AML_DAY1_N2)), rep("B.Day7", ncol(AML_DAY7_N1)), rep("B.Day7", ncol(AML_DAY7_N2)), rep("C.Day14", ncol(AML_DAY14_N1)), rep("C.Day14", ncol(AML_DAY14_N2)))
# 3- Use Asaf's function to sum up duplicate rows - SAVE FILE AFTER

dups=fData(AML_DAY1_N1)$symbol[which(duplicated(fData(AML_DAY1_N1)$symbol))]; rmv=NULL
for (i in dups)
{
  ind=which(fData(AML_DAY1_N1)$symbol==i)
  minInd=min(ind); ind=setdiff(ind,minInd)
  if(length(ind)==1) all[minInd,]=all[minInd,]+all[ind,]
  else all[minInd,]=all[minInd,]+colSums(all[ind,])
  rmv=c(rmv,ind)
  print(match(i,dups))
}
all=all[-rmv,]
rownames(all)=fData(AML_DAY1_N1)$symbol[match(rownames(all),fData(AML_DAY1_N1)$id)]

save(all, file = "~/Desktop/AML_Eirini_Josi/fullUmitab.Rdata")

source("~/Desktop/Single cell project/R studio analysis/LogNorm.R")
all=all[rowSums(all)>0,]
data=all
#data=as.matrix(all[,sampleType %in% c("wt1","wt2")]) # if you want to only look at a subset
#data=data[rowSums(data)>0,] #not all genes are expressed, remove empty rows

##Change variable gene detection to effi's approach
variableGenes=function(dataMat,xmincutoff=NULL,ymincutoff=NULL,dsTo=NULL,seedNum=10021){
  set.seed(seedNum)
  if(is.null(dsTo)) dsTo=min(colSums(dataMat))
  downsamp_one=function(x,n){
    hist(sample(rep(1:length(x),x),size = n,replace=F),breaks=(0:length(x)+.5),plot=F)$counts}
  a=dataMat[,colSums(dataMat)>=dsTo];message(dim(a)[2])
  ds=apply(a,2,downsamp_one,dsTo)
  rownames(ds)=rownames(dataMat)
  ds=ds[rowMeans(ds)>0,]
  means=rowMeans(ds)
  if(is.null(xmincutoff)) xmincutoff=quantile(log(means))[3]
  varGenes=apply(ds,1,var)/means
  if(is.null(ymincutoff)) ymincutoff=quantile(log(varGenes))[3]
  plot(log(means),log(varGenes))
  abline(h = ymincutoff)
  abline(v = xmincutoff)
  varGenesNames=rownames(ds)[which(log(varGenes)>ymincutoff & means >xmincutoff )]
  return(varGenesNames)
}

varGenes = variableGenes(dataMat = all, xmincutoff = NULL, ymincutoff = 0.25, dsTo = 5000, seedNum = 100)
save(varGenes, file = "~/Desktop/AML_Eirini_Josi/varGenes_ds500_xminNULL_ymin025.Rdata")
barcodes = colnames(all)
save(barcodes, file = "~/Desktop/AML_Eirini_Josi/barcodes")

cellNames = c(paste0("DAY1_N1_",1:ncol(AML_DAY1_N1)), paste0("DAY1_N2_",1:ncol(AML_DAY1_N2)), paste0("DAY7_N1_",1:ncol(AML_DAY7_N1)), paste0("DAY7_N2_",1:ncol(AML_DAY7_N2)), paste0("DAY14_N1_",1:ncol(AML_DAY14_N1)), paste0("DAY14_N2_",1:ncol(AML_DAY14_N2)))
colnames(all) = cellNames
names(mapping) = colnames(all)
names(day) = colnames(all)
save(mapping, file = "~/Desktop/AML_Eirini_Josi/mapping.Rdata")
save(day, file = "~/Desktop/AML_Eirini_Josi/day.Rdata")
