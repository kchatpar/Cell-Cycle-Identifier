#!/usr/bin/env Rscript

#Merges the UMItabs and splits by each cluster.
library(ggplot2)
library(matrixStats)
setwd("/gpfs/commons/home/chatpark-1287/Cell_Cycle_Plot2/")

#Loads the data.
load("/Users/krishnachatpar/Downloads/Franco10X_all.RData")
#load("~/Cell_Cycle_Data/Franco10X_all.RData")

#Split UMItab into different genotypes
wt1_mat <- as.matrix(all[ , grepl( "wt1" , colnames( all ) ) ])
wt2_mat <- as.matrix(all[ , grepl( "wt2" , colnames( all ) ) ])
Idh2_mat <- as.matrix(all[ , grepl( "Idh2" , colnames( all ) ) ])
Dnmt_mat <- as.matrix(all[ , grepl( "Dnmt3a2" , colnames( all ) ) ])
print("Load Succesful")

#Filter Tab's
wt1_mat_filt = wt1_mat[rowSums(wt1_mat)>0,]
dim(wt1_mat_filt)
wt2_mat = wt2_mat[rowSums(wt2_mat)>0,]
Idh2_mat = Idh2_mat[rowSums(Idh2_mat)>0,]
Dnmt_mat = Dnmt_mat[rowSums(Dnmt_mat)>0,]

#Normalize matricies
norm_wt1<-apply(X = wt1_mat_filt[,],MARGIN = 2,FUN = function(x) x/sum(x))

norm_wt2<-apply(X = wt2_mat[,],MARGIN = 2,FUN = function(x) x/sum(x))
norm_Idh2<-apply(X = Idh2_mat[,],MARGIN = 2,FUN = function(x) x/sum(x))
norm_Dnmt<-apply(X = Dnmt_mat[,],MARGIN = 2,FUN = function(x) x/sum(x))

norm_wt1_df<-as.data.frame(norm_wt1)
norm_wt2_df<-as.data.frame(norm_wt2)
norm_Idh2_df<-as.data.frame(norm_Idh2)
norm_Dnmt<-as.data.frame(norm_Dnmt)

split_wt<-list(norm_wt1_df,norm_wt2_df)
split_ko<-list(norm_Idh2,norm_Dnmt)

names(split_wt) = c("wt_1", "wt_2")
names(split_ko) = c("Idh2", "Dnmt")

master.G1 = list()
master.G2 = list()

c1_df<-list()

#Loop through wt matricies and create updated gene lists.
for(i in 1:length(split_wt)){
  #Load G1_S genes
  Input_G1_S_Genes<-read.delim("/Users/krishnachatpar/Desktop/New York Genome Center:Cornell/Cell Cycle Data/G1_S_Genes.txt",header = F,sep=" ")
  G1<-t(Input_G1_S_Genes)
  G1_S_Genes<-as.vector(G1[,1])

  #Load the G2_M genes and store it in G2_M_Genes
  Input_G2_M_Genes<-read.delim("/Users/krishnachatpar/Desktop/New York Genome Center:Cornell/Cell Cycle Data/G2_M_Genes.txt",header = F,sep=" ")
  G2<-t(Input_G2_M_Genes)
  G2_M_Genes<-as.vector(G2[,1])

  #Go through wt1 and wt2 and make a master gene list
  c1_wt = split_wt$wt_1
  dim(c1_wt)
  
  #Transpose and Correlate Genes
  c1_t_wt<-t(c1_wt)
  c1_df[[i]]<-as.data.frame(cor(c1_t_wt))

  #Starts the updating of the gene lists
  new.genes = list()
  #Find genes with correlation > 0.90 for G1_S Phase
  c1_df_filtered = c1_df[[i]][intersect(G1_S_Genes, rownames(c1_df[[i]])),]
  
  #New Genes2 = correlated genes above a certain threshold.
  new.genes = apply(c1_df_filtered, MARGIN = 1, function(x) x = x[which(x > 0.90)])
  new.genes2 = lapply(new.genes, function(x) names(x))
  
  #New Gene List for G1_S = genes.updated
  genes.updated = unique(c(unlist(new.genes2), G1_S_Genes))
  master.G1[[i]] = genes.updated
  
  #Find genes with cor value>0.90 for G2_M phase
  new.genes2 = list()
  c1_df_filtered_G2 = c1_df[[i]][intersect(G2_M_Genes, rownames(c1_df[[i]])),]
  new.genes = apply(c1_df_filtered_G2, MARGIN = 1, function(x) x = x[which(x > 0.90)])
  new.genes2 = lapply(new.genes, function(x) names(x))
  
  #New Gene List G2_M = genes.updated_G2
  genes.updated_G2 = unique(c(unlist(new.genes2), G2_M_Genes))
  master.G2[[i]] = genes.updated_G2
}

G1_Intersect<-master.G1[[1]]
G2_Intersect<-master.G2[[1]]

for(i in 1:length(master.G1)){
  G1_Intersect<-union(G1_Intersect,master.G1[[i]])
  G2_Intersect<-union(G2_Intersect,master.G2[[i]])
}

#Updates G1_S Gene List
G1_S_Genes<-G1_Intersect

#Updates G2_M Genes
G2_M_Genes<-G2_Intersect

#Remove unique genes from G1/S, G2/M Lists

temp_G1_1<-G1_S_Genes
temp_G1_2<-G1_S_Genes

temp_G2_1<-G2_M_Genes
temp_G2_2<-G2_M_Genes

temp_G1_diff<-setdiff(temp_G1_1,temp_G2_1)
temp_G2_diff<-setdiff(temp_G2_2,temp_G1_2)

G1_S_Genes<-temp_G1_diff
G2_M_Genes<-temp_G2_diff

#Generate mean names for files
wt1_file_name<-"Wt1_vs_Idh2_Mean.png"

#plot_mean_name<-c()
#for(i in 1:length(split_wt)){
#  plot_mean_name[i]<-paste(clustAnnots[[i]],"mean.png",sep="_")
#}  

#Generate sum names for files
wt1_file_name<-"Wt1_vs_Idh2_Sum.png"

#plot_sum_name<-c()
#for(i in 1:length(hsc_clust)){
#  clust<-clustAnnots[[i]]
#  plot_sum_name[i]<-paste(clustAnnots[[i]],"sum.png",sep="_")
#}

#for(i in 1:length(hsc_clust)){
#Split based on type

#Load for Idh2 Umitab
  c1_tet = c1_tet[rowSums(c1_tet)>0,]
  norm_c1<-split_ko$Idh2

#Load cluster 1 for wt
  norm_c1_wt<-split_wt$wt1

#Normalize ko vs. wt
   norm_ko_df<-split_ko$Idh2
   norm_wt_df<-split_wt$wt1

#Find intersecting genes
  G1_S_ko_df<-norm_ko_df[intersect(rownames(norm_ko_df),G1_S_Genes),]
  G1_S_wt_df<-norm_wt_df[intersect(rownames(norm_wt_df),G1_S_Genes),]

  G2_M_ko_df<-norm_ko_df[intersect(rownames(norm_ko_df),G2_M_Genes),]
  G2_M_wt_df<-norm_wt_df[intersect(rownames(norm_wt_df),G2_M_Genes),]

#Sum each cell.
  G1_S_ko_sum_df<-colSums(G1_S_ko_df)
  G1_S_ko_sum_df<-G1_S_ko_sum_df/max(G1_S_ko_sum_df)

  G1_S_wt_sum_df<-colSums(G1_S_wt_df)
  G1_S_wt_sum_df<-G1_S_wt_sum_df/max(G1_S_wt_sum_df)

  G2_M_ko_sum_df<-colSums(G2_M_ko_df)
  G2_M_ko_sum_df<-G2_M_ko_sum_df/max(G2_M_ko_sum_df)

  G2_M_wt_sum_df<-colSums(G2_M_wt_df)
  G2_M_wt_sum_df<-G2_M_wt_sum_df/max(G2_M_wt_sum_df)

#Mean of each cell
   #G1_S_ko_mean_df<-colMeans(G1_S_ko_df)
   #G1_S_ko_mean_df<-G1_S_ko_mean_df/max(G1_S_ko_mean_df)

   #G1_S_wt_mean_df<-colMeans(G1_S_wt_df)
   #G1_S_wt_mean_df<-G1_S_wt_mean_df/max(G1_S_wt_mean_df)

   #G2_M_ko_mean_df<-colMeans(G2_M_ko_df)
   #G2_M_ko_mean_df<-G2_M_ko_mean_df/max(G2_M_ko_mean_df)

   #G2_M_wt_mean_df<-colMeans(G2_M_wt_df)
   #G2_M_wt_mean_df<-G2_M_wt_mean_df/max(G2_M_wt_mean_df)

#Create DF's for Ko vs. Wt for the sum
   df1 = as.data.frame(cbind(G1_S_wt_sum_df, G2_M_wt_sum_df))
   colnames(df1) = c("G1_S", "G2_M")

   df2 = as.data.frame(cbind(G1_S_ko_sum_df, G2_M_ko_sum_df))
   colnames(df2) = c("G1_S", "G2_M")

   df = as.data.frame(rbind(df1,df2))
   df$Genotype = c(rep("Wt1", length(G1_S_wt_sum_df)), rep("Idh2 KO", length(G1_S_ko_sum_df)))

#Create Df's for Ko vs. Wt for the mean
   #df1_mean = as.data.frame(cbind(G1_S_wt_mean_df, G2_M_wt_mean_df))
   #colnames(df1_mean) = c("G1_S", "G2_M")

   #df2_mean = as.data.frame(cbind(G1_S_ko_mean_df, G2_M_ko_mean_df))
   #colnames(df2_mean) = c("G1_S", "G2_M")

   #df_mean = as.data.frame(rbind(df1_mean,df2_mean))
   #df_mean$Genotype = c(rep("Wt", length(G1_S_wt_mean_df)), rep("Tet2 KO", length(G1_S_ko_mean_df)))

   #p = ggplot(data = df_mean, aes(x = G1_S, y = G2_M, colour = Genotype))+
    # geom_point(alpha = 0.5)+
    # theme_classic()+
    # labs(title = paste0(clustAnnots[n]), y = "G2_M mean", x = "G1_S mean")
   #print(p)

  #ggsave()

   q = ggplot(data = df, aes(x = G1_S, y = G2_M, colour = Genotype))+
     geom_point(alpha = 0.5)+
     theme_classic()+
     labs(title = paste0(clustAnnots[n]), y = "G2_M sum", x = "G1_S sum")
   print(q)
   ggsave(wt1_file_name)
  #}
  
  #***Heat Map Plotting of Correlation***
  #  heatmap(as.matrix(c1_df_filtered), scale="column",col=heat.colors(256),main=plot_hm_name[i],Rowv = NA,Colv = NA)
  #  as.matrix(c1_df_filtered)
  #  dim(c1_df_filtered)
  #  png(plot_hm_name[i])
  #  heatmap(as.matrix(c1_df_filtered), scale="column",col=heat.colors(256),main=clustAnnots[[n]],Rowv = NA,Colv = NA)
  #  dev.off()
  #}

  #Create file names for plot saving
  plot_hm_name<-c()
  for(i in 1:length(split_all)){
    clust<-split_all[[i]]
    plot_hm_name[i]<-paste(split_all[[i]],"heat_map.png",sep="_")
  }



