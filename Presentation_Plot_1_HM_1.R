#Convert human list to mouse list to test on Nir's data.
library(ggplot2)
library(matrixStats)
library(Matrix)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(grid)

#Load Control Data
load("~/Desktop/New York Genome Center:Cornell/Cell Cycle Data/nir10x_monoMacroTP2.RData")

#Load the M genes and store it in M_Genes
Input_M_Genes<-read.delim("~/Desktop/New York Genome Center:Cornell/10X_Analysis/10X_Gene_Lists/Human_M_Gene_List.txt",header = F,sep=",")
M_t<-t(Input_M_Genes)
M_Genes<-as.vector(M_t[,1])

#Load S Genes
Input_S_Genes<-read.delim("~/Desktop/New York Genome Center:Cornell/10X_Analysis/10X_Gene_Lists/Human_S_Gene_List.txt",header = F,sep=",")
S_t<-t(Input_S_Genes)
S_Genes<-as.vector(S_t[,1])

#Load G1 Genes ***New Gene List***
Input_G1_Genes<-read.delim("~/Desktop/New York Genome Center:Cornell/10X_Analysis/10X_Gene_Lists/Human_G1_Gene_List.txt",header = F,sep=",")
G1_t<-t(Input_G1_Genes)
G1_Genes<-as.vector(G1_t[,1])

#Load G2 Genes
Input_G2<-read.delim("~/Desktop/New York Genome Center:Cornell/10X_Analysis/10X_Gene_Lists/Human_G2_Gene_List.txt",header = F,sep=",")
G2_t<-t(Input_G2)
G2_Genes<-as.vector(G2_t[,1])

#Start conversion of genes.
G1_Genes<-convertHumanGeneList(G1_Genes)
G2_Genes<-convertHumanGeneList(G2_Genes)
M_Genes<-convertHumanGeneList(M_Genes)
S_Genes<-convertHumanGeneList(S_Genes)
length(G1_Genes)
length(G2_Genes)
length(M_Genes)
length(S_Genes)

#Split UMItab into different genotypes
wt1_mat <- data

#Filter Tab's
wt1_mat_filt = wt1_mat[rowSums(wt1_mat)>0,]

G1_S_Genes<-G1_Genes
G2_M_Genes<-G2_Genes

G1_S_Genes_mat<-wt1_mat[intersect(G1_S_Genes, rownames(wt1_mat)),]
g1_rsum<-c(rowSums(G1_S_Genes_mat))
filt_g1<-g1_rsum[g1_rsum<2000]
G1_S_Genes<-names(filt_g1)
length(G1_S_Genes)

G2_df<-wt1_mat[intersect(G2_M_Genes, rownames(wt1_mat)),]
g2_rsum<-c(rowSums(G2_df))
filt_g2<-g2_rsum[g2_rsum<2000]
G2_M_Genes<-names(filt_g2)
length(G2_M_Genes)

S_df<-wt1_mat[intersect(S_Genes, rownames(wt1_mat)),]
S_rsum<-c(rowSums(S_df))
filt_S<-S_rsum[S_rsum<2000]
S_Genes<-names(filt_S)

M_df<-wt1_mat[intersect(M_Genes, rownames(wt1_mat)),]
M_rsum<-c(rowSums(M_df))
filt_M<-M_rsum[M_rsum<2000]
M_Genes<-names(filt_M)
length(M_Genes)

#Normalize matricies
norm_wt1<-apply(X = wt1_mat_filt,MARGIN = 2,FUN = function(x) x/sum(x))

#Subset norm1
norm_wt1<-norm_wt1[,sample(x = c(1:ncol(norm_wt1)), size = 1000, replace = F)]

#Get Unique genes minus the G0 phase.
geneNames = c(G1_S_Genes,G2_M_Genes,S_Genes,M_Genes)
geneFactor = c(rep("G1",length(G1_S_Genes)),rep("G2",length(G2_M_Genes)), rep("M",length(M_Genes)),rep("S",length(S_Genes)))
labels = c("G1","G2","M","S")
names(geneFactor) = geneNames
uniqueGenes = unique(geneNames)
geneFactor = geneFactor[uniqueGenes]

annotation_df = as.data.frame(cbind(uniqueGenes,geneFactor))
rownames(annotation_df) = annotation_df$uniqueGenes
annotation_df$uniqueGenes = NULL

#Subset normalized data frame for unique genes and filter genes with no data.
#Remove Mitochondiral Genes.

unique_norm_df = norm_wt1[uniqueGenes,]
mito.genes <- grep("^MT-", rownames(norm_wt1), value = T, ignore.case = T)
unique_norm_df<-unique_norm_df[!rownames(unique_norm_df) %in% mito.genes,]
unique_norm_df = unique_norm_df[rowSums(unique_norm_df)>0,]
unique_norm_df = unique_norm_df[rowSds(unique_norm_df)>.00005,]

pheatmap(
  mat    = log10(unique_norm_df+1e-5),
  color  = inferno(15),
  #breaks = mat_breaks,
  #border_color = NA,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_row = annotation_df,
  cluster_rows = T,
  scale = "row"
  #display_numbers = T
  #drop_levels = TRUE,
  #fontsize = 14
  #kmeans_k = 4
)
