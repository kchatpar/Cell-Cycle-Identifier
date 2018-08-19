library("Matrix")
library("monocle")
library("Biobase")
library("stringi")
library("reshape2")
library(plyr)
library(igraph)
library(beepr)

setwd(dir = "~/Desktop/AML_Eirini_Josi/")
#load("mapping.Rdata")
#load("varGenes_ds500_xminNULL_ymin025.Rdata")
#load("aml_Seurat_object.Rdata")
#load("allMarkers.Rdata")
load("top_markers.Rdata")
set.seed(seed = 16849)

# 1 - Randomly sample cells from original data

random.vector = sample(x = c(1:ncol(aml@raw.data)), size = 5000, replace = F)
cells = aml@raw.data[,random.vector]

# 2 - Generate annotations and monocle object containing the raw data

sampleType = as.data.frame(mapping[colnames(cells)])
colnames(sampleType) = "Condition"
pd=new("AnnotatedDataFrame", data = sampleType)
                     
AML <- newCellDataSet(as.matrix(cells), phenoData = pd, expressionFamily = negbinomial())

# 3 - Estimate variable genes within data
AML <- estimateSizeFactors(AML)
AML <- estimateDispersions(AML)
print(head(pData(AML)))

pData(AML)$Total_mRNAs <- Matrix::colSums(exprs(AML))
AML <- AML[,pData(AML)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(AML)$Total_mRNAs)) + 2*sd(log10(pData(AML)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(AML)$Total_mRNAs)) - 3*sd(log10(pData(AML)$Total_mRNAs)))
qplot(Total_mRNAs, data=pData(AML), geom="density") +
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)

# 4 - Select variable genes

AML <- AML[,pData(AML)$Total_mRNAs > lower_bound &
               pData(AML)$Total_mRNAs < upper_bound]
AML <- detectGenes(AML, min_expr = 0.1)
disp_table <- dispersionTable(AML)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.02 & dispersion_empirical >= 1.25 * dispersion_fit)
AML <- setOrderingFilter(AML, unsup_clustering_genes$gene_id)
plot_ordering_genes(AML)
AML <- clusterCells(AML, num_clusters= 7)

# 5 - Substract phenodata from Seurat analysis and keep expressed genes

pData(AML)$Cluster=aml@data.info[rownames(pData(AML)),]$cluster_name
pData(AML)$Day = aml@data.info[rownames(pData(AML)),]$day

# 6 - Get marker genes for semi-supervised pseudotime ordering

AML_expressed_genes <- row.names(subset(fData(AML), num_cells_expressed >= 5))
AML_filtered <- AML[AML_expressed_genes,]
exprs_filtered <- t(t(exprs(AML_filtered)/pData(AML_filtered)$Size_Factor))
nz_genes <- which(exprs_filtered != 0)
exprs_filtered[nz_genes] <- log(exprs_filtered[nz_genes] + 1)
# 6 - Calculate the variance across genes without converting to a dense matrix:
expression_means <- Matrix::rowMeans(exprs_filtered)
expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
# 7 - Filter out genes that are constant across all cells:
genes_to_keep <- expression_vars > 0
exprs_filtered <- exprs_filtered[genes_to_keep,]
expression_means <- expression_means[genes_to_keep]
expression_vars <- expression_vars[genes_to_keep]
# Here's how to take the top PCA loading genes, but using
# sparseMatrix operations the whole time, using irlba.
irlba_pca_res <- irlba(t(exprs_filtered),
                       nu=0,
                       center=expression_means,
                       scale=sqrt(expression_vars),
                       right_only=TRUE)$v
row.names(irlba_pca_res) <- row.names(exprs_filtered)
PC1_genes <- names(sort(abs(irlba_pca_res[, 1]), decreasing = T))[1:50]
PC2_genes <- names(sort(abs(irlba_pca_res[, 2]), decreasing = T))[1:50]
PC3_genes <- names(sort(abs(irlba_pca_res[, 3]), decreasing = T))[1:50]
ordering_genes <- union(union(PC1_genes, PC2_genes),PC3_genes)

#If we want to remove Gm and ribosomal genes
ordering_genes <- ordering_genes[grep(ordering_genes, pattern = "Gm", invert = T)]
ordering_genes <- ordering_genes[grep(ordering_genes, pattern = "RPS", invert = T)]
ordering_genes <- ordering_genes[grep(ordering_genes, pattern = "RPL", invert = T)]
ordering_genes <- ordering_genes[grep(ordering_genes, pattern = "RP", invert = T)]

AML <- setOrderingFilter(AML, top_markers$gene)
AML <- reduceDimension(AML, max_components=2)
AML <- orderCells(AML, reverse=FALSE)

plot_cell_trajectory(AML, color_by= "Day",show_branch_points = F, show_tree = T, show_backbone = T, cell_size = 0.25)+
  scale_color_manual(values = c("#E53F0D", "#6DBE44", "#4E81C2"))+
  facet_wrap(~Day)
 
plot_cell_trajectory(AML, color_by= "Cluster",show_branch_points = T, show_tree = T, show_backbone = T, cell_size = 0.25)

plot_cell_trajectory(AML, color_by= "Cluster",show_branch_points = F, show_tree = T, show_backbone = T, cell_size = 0.25)+
  facet_wrap(~Day)

plot_cell_trajectory(AML, color_by = "State", show_branch_points = F, cell_size = 0.25) + facet_wrap(~State)
AML <- orderCells(AML, root_state= 34) # Correct initial state

plot_cell_trajectory(AML, color_by = "Pseudotime", show_branch_points = F, cell_size = 0.25) + 
  scale_color_continuous(low = c("#F9A71B"), high = c("#F2111C")) +
  theme_classic(base_size = 16)


plot_cell_trajectory(HSMM, color_by = "CellType", show_branch_points = F) + facet_wrap(~CellType) +
  scale_color_manual(values = c("#F9A61B","#F2111C", "#F2111C", "#7094CB", "#7094CB", "#7094CB", "#F9A61B", "#F9A61B", "#7094CB"))

plot_cell_trajectory(HSMM, color_by = "Origin", show_branch_points = F, show_tree = F,cell_size = 0.25 ) + facet_wrap(~Origin ) +
  scale_color_manual(values = c("#7094CB","#F3756B"))

plot_cell_trajectory(HSMM, color_by = "Origin", show_branch_points = F, show_tree = F,cell_size = 0.25 ) +
  scale_color_manual(values = c("#7094CB","#F3756B"))


plot_cell_trajectory(HSMM, color_by = "CellType", show_branch_points = F, show_tree = F) +
  scale_color_manual(values = c("#F9A61B","#F2111C", "#F2111C", "#7094CB", "#7094CB", "#7094CB", "#F9A61B", "#F9A61B", "#7094CB"))+
  theme_classic(base_size = 16)

plot_genes_branched_pseudotime(cds = AML, branch_point = 10, color_by = "Cluster", ncol = 2) +
  theme_classic(base_size = 16)

# - To plot the expression of adhesion molecules across trajectories

Genes_from_bulk <- as.data.frame(read_csv("~/Desktop/AML_Eirini_Josi/Gene_lists/Genes_from_bulk.csv", col_names = FALSE))
Genes_from_bulk = Genes_from_bulk$X1

plot_genes_branched_pseudotime(cds = AML[intersect(Genes_from_bulk, rownames(AML)),], branch_point = 10, color_by = "Day", ncol = 4)+
  theme_classic(base_size = 16)


BEAM_res = BEAM(HSMM, branch_point = 10, cores = 1)
BEAM_res = BEAM_res[order(BEAM_res$qval),]
plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 10,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)



