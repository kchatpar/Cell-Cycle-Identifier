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
library(cellrangerRkit)
library(beepr)
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R") - To install cellRanger kit
library(cellrangerRkit)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(readr)

setwd(dir = "~/Desktop/AML_Eirini_Josi/")
set.seed(seed = 4123)

load("fullUmitab.Rdata")
load("mapping.Rdata")
load("varGenes_ds500_xminNULL_ymin025.Rdata")
load("aml_Seurat_object.Rdata")
load("allMarkers.Rdata")
# 1- QC metrics

umis = as.data.frame(cbind(colSums(all), mapping))
colnames(umis) = c("umis", "mapping")
umis$umis = as.numeric(as.character(umis$umis))

detectedGenes = apply(all, MARGIN = 2, function(x) length(which(x>0)))
genes = as.data.frame(cbind(detectedGenes, mapping))
colnames(genes) = c("genes", "mapping")
genes$genes = as.numeric(as.character(genes$genes))

repDay1 = as.data.frame(cbind(rowSums(AML_DAY1_N1), rowSums(AML_DAY1_N2)))
repDay1 = log(repDay1, base = 2)
cor(repDay1$V1, repDay1$V2, method = "spearman")
ggplot(data = repDay1, aes(x = V1, y = V2))+
  geom_point(color = "#E53F0C", alpha = 0.5)+theme_classic(base_size = 16)+
  labs(x = "Replicate 1", y = "Replicate 2", title = "Day 1")

repDay7 = as.data.frame(cbind(rowSums(AML_DAY7_N1), rowSums(AML_DAY7_N2)))
repDay7 = log(repDay7, base = 2)
cor(repDay7$V1, repDay7$V2, method = "spearman")
ggplot(data = repDay7, aes(x = V1, y = V2))+
  geom_point(color = "#6DBE45", alpha = 0.5)+theme_classic(base_size = 16)+
  labs(x = "Replicate 1", y = "Replicate 2", title = "Day 7")

repDay14 = as.data.frame(cbind(rowSums(AML_DAY14_N1), rowSums(AML_DAY14_N2)))
repDay14 = log(repDay14, base = 2)
cor(repDay14$V1, repDay14$V2, method = "spearman")
ggplot(data = repDay14, aes(x = V1, y = V2))+
  geom_point(color = "#4E81C2", alpha = 0.5)+theme_classic(base_size = 16)+
  labs(x = "Replicate 1", y = "Replicate 2", title = "Day 14")

ggplot(data = umis, aes(x = mapping, y = umis, fill = mapping)) + 
  geom_jitter(alpha = 0.05, width = 0.1, aes(color = mapping))+
  geom_boxplot(outlier.alpha = 0, notch = T) + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = c("#E53F0C", "#E86036", "#6DBE45", "#A6E085", "#4E81C2", "#7FAFE2"))+
  scale_color_manual(values = c("#E53F0C", "#E86036", "#6DBE45", "#A6E085", "#4E81C2", "#7FAFE2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Numer of UMIs per cell", x = "")

ggplot(data = genes, aes(x = mapping, y = genes, fill = mapping)) + 
  geom_jitter(alpha = 0.05, width = 0.1, aes(color = mapping))+
  geom_boxplot(outlier.alpha = 0, notch = T) + 
  theme_classic(base_size = 16) + 
  scale_fill_manual(values = c("#E53F0C", "#E86036", "#6DBE45", "#A6E085", "#4E81C2", "#7FAFE2"))+
  scale_color_manual(values = c("#E53F0C", "#E86036", "#6DBE45", "#A6E085", "#4E81C2", "#7FAFE2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Numer of detected genes per cell", x = "")

# 2- Run Seurat for clustering and visualization

aml <- new("seurat", raw.data =as.matrix(all))
aml <- Setup(aml, min.cells = 3, min.genes = 200, do.logNormalize = F, total.expr = 1e5, project = "AML_Eirini_Josi")
aml@data=LogNormalize(aml@data,scale.factor = 1e5)

# 3- Perform dimension reduction by PCA
aml <- PCA(aml, pc.genes = varGenes, do.print = TRUE, pcs.print = 5, genes.print = 5)
aml <- ProjectPCA(aml)

# 4- Decide the number of PCs for downstream analysis
PCElbowPlot(aml,num.pc = 40)
PCHeatmap(aml, pc.use =1:10 , cells.use = 100, do.balanced = TRUE)

aml <- FindClusters(aml, pc.use = 1:20, resolution = 0.25, print.output = 0, save.SNN = T)

aml <- RunTSNE(object = aml, dims.use = 1:40, do.fast = T, dim_embed = 2, reduction.use = "pca")
TSNEPlot(aml, do.label = T, pt.size = 1)

#Find cluster markers
amlMarkers <- FindAllMarkers(aml, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.3)
save(amlMarkers, file = "allMarkers.Rdata")
D1vsD7.cell.markers <- FindMarkers(aml,ident.1 = "DAY1",ident.2 = "DAY7",only.pos = FALSE, min.pct = 0.25, thresh.use = 0.4)
D7vsD14.cell.markers <- FindMarkers(aml,ident.1 = "DAY1",ident.2 = "DAY7",only.pos = FALSE, min.pct = 0.25, thresh.use = 0.4)

save(aml, file = "aml_Seurat_object.Rdata")

top_markers <- amlMarkers %>% group_by(cluster) %>% top_n(15, avg_diff)
write_csv(x = top_markers, path = "Top50_genes_r025.csv")
View(top_markers)
save(top_markers, file = "top_markers.Rdata")


heatmapData = as.matrix(aml@data[intersect(top_markers$gene, rownames(aml@data)),sample(x = c(1:ncol(aml@data)), size = 3000, replace = F)])

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
clustAnnotation = as.data.frame(aml@data.info$res.0.25)
mat_breaks <- quantile_breaks(heatmapData, n = 11)
pheatmap(
  mat               = heatmapData,
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Quantile Color Scale"
)

tsneCoords = aml@tsne.rot
tsneCoords$Condition = mapping[rownames(tsneCoords)]
tsneCoords$Day = day[rownames(tsneCoords)]
tsneCoords$Cluster = aml@data.info$res.0.25

ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = Cluster)) +
  geom_point(size = 0.5) 

ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = Cluster)) +
  geom_point(size = 0.25) +
  facet_wrap(~Day, ncol = 1) 


ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = Cluster)) +
  geom_point(size = 0.5) +
  facet_wrap(~Condition, ncol = 2) 
tsneCoords$MPO = aml@data["MPO",]

gene = top_markers$gene
for(i in 1:length(gene)){
p = ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = aml@data[gene[i],])) +
  geom_point(size = 0.5) +
  scale_color_continuous(low = "#FEDE00", high = "#E53F0C")+
  labs(color = paste0(gene[i], " expression"), title = paste0("Cluster ",top_markers$cluster[i]," marker gene"))
print(p)
}

# Cells in cluster by day
cellNumbers = as.data.frame.matrix(table(tsneCoords$Condition, tsneCoords$Cluster))
cellFrequencies = apply(cellNumbers, MARGIN = 1, function(x) x/sum(x))

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(as.matrix(cellFrequencies), n = 11)
pheatmap(
  mat               = as.matrix(cellFrequencies),
  color             = inferno(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Quantile Color Scale"
)

# Batch
batch = cor(cellFrequencies)

plotData = melt(cellFrequencies)

ggplot(PlotData, aes(x= Var2, y=MeanFrequency, group = Genotype)) + 
  geom_errorbar(aes(ymin=MeanFrequency-SdFrequency, ymax=MeanFrequency+SdFrequency), width=.1, position = position_dodge(0.9)) +
  geom_bar(position = position_dodge(), stat = "identity", aes(fill = Genotype))+
  scale_fill_manual(values = c("#7194CA","#F3756B"))+
  labs(y = "Mean frequency", x = "Cluster")+
  theme_classic(base_size = 16) +
  ylim(0,0.12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Cluster 0: MEG/HSC
#Cluster 1: Mast cell/HPC
#Cluster 2: stromal cells
#Cluster 3: Lymphocyte B+T
#Cluster 4: Down in LSCs
#Cluster 5: Up in HSPCs
#Cluster 6: hESCs

current.cluster.ids = c("0", "1", "2", "3", "4", "5", "6")
new.cluster.ids = c("MEG/HSC", "Mst/HPC", "Stromal cells", "Lymphocyte B+T", "Down in LSCs", "Up in HSPCs", "hESCs")
aml@ident <- plyr::mapvalues(aml@ident, from = current.cluster.ids, to = new.cluster.ids)
aml@data.info$cluster_name = aml@ident
aml@data.info$day = rep("undet", nrow(aml@data.info))
aml@data.info$day[grep(rownames(aml@data.info), pattern = "DAY1")] = "A-DAY1"
aml@data.info$day[grep(rownames(aml@data.info), pattern = "DAY7")] = "B-DAY7"
aml@data.info$day[grep(rownames(aml@data.info), pattern = "DAY14")] = "C-DAY14"

save(aml, file = "aml_Seurat_object.Rdata")

# - To generate scores for expression of adhesion molecules

Adhesion_gene_list <- as.data.frame(read_csv("~/Desktop/AML_Eirini_Josi/Gene_lists/Adhesion_gene_list.csv"))
Adhesion_gene_list = Adhesion_gene_list$Adhesion_genes

stickyScore = colSums(aml@data[intersect(Adhesion_gene_list, rownames(aml@data)),])

tsneCoords$AdhesionScore = stickyScore

ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = AdhesionScore)) +
  geom_point(size = 0.5) +
  scale_color_continuous(low = "#EDEDF2", high = "#000CE2")+
  facet_wrap(~aml@data.info$day)

plot(hist(rowSums(aml@data[intersect(Adhesion_gene_list, rownames(aml@data)),])))

Adhesion = aml@data[intersect(Adhesion_gene_list, rownames(aml@data)),]
Adhesion = Adhesion[rowSums(Adhesion) > 10000,]
highExpressedAdhesion = rownames(Adhesion)

for(i in 1:length(highExpressedAdhesion)){
p = ggplot(data = tsneCoords, aes( x = tSNE_1, y = tSNE_2, color = aml@data[highExpressedAdhesion[i],])) +
  geom_point(size = 0.5) +
  scale_color_continuous(low = "#FEDE00", high = "#E53F0C")+
  facet_wrap(~aml@data.info$day)+
  labs(title = paste0(highExpressedAdhesion[i], " expression"))
print(p)
}









