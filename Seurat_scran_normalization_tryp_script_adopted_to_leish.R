library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(MAST)
library(scater)
library(scran)
library(clustree)
library(RColorBrewer)

# Load data set
tryp.data <- Read10X(data.dir = "D:/GCRF_UoG/Felix_raw_data/leish_dset/")

# Filter out VSG (need to import text file VSG_list - VSGs from TriTryp search)
# VSGs <-as.character(VSG_list$V1)
# tryp.data <- tryp.data[! rownames(tryp.data) %in% VSGs,]

# Initialize the Seurat object with the raw (non-normalized data).
tryp_raw_object <- CreateSeuratObject(counts = tryp.data, project = "Cells", min.cells = 5, min.features = 3)

##Cell death genes
mt_marker_genes <- c("ND7", "COIII", "MURF4", "MURF1", "ND1", "MURF2", "COI", "COII", "ND4", "RPS12", "ND5", "Cyb")
tryp_raw_object[["percent.MT"]] <- PercentageFeatureSet(tryp_raw_object, pattern = "^MT-")

# Visualize QC metrics as a violin plot, need way of using mito genes
featureVln <- VlnPlot(tryp_raw_object, features = "nFeature_RNA", pt.size = 0.01) + labs(title = "Features") + ylab("Count") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
UMIVln <- VlnPlot(tryp_raw_object, features = "nCount_RNA", pt.size = 0.01) + labs(title = "UMIs") + ylab("Count") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
MTVln <- VlnPlot(tryp_raw_object, features = "percent.MT", pt.size = 0.01) + labs(title = "kDNA Genes") + ylab("% Features") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 10), axis.text.y =  element_text(size = 10))
ggarrange(featureVln, UMIVln, MTVln,ncol = 3, legend = "none", label.x = 0)

UMI_MT <- FeatureScatter(tryp_raw_object, feature1 = "nCount_RNA", feature2 = "percent.MT", pt.size = 0.1, smooth = FALSE) + labs(title = "") + scale_y_continuous(name = "% kDNA genes") + xlab("UMI Count") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
UMI_feature <- FeatureScatter(tryp_raw_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)+ labs(title = "") + ylab("Feature Count") + xlab("UMI Count") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
ggarrange(UMI_MT, UMI_feature,  ncol = 2, legend = "none")

ggplot(tryp_raw_object@meta.data, aes(x = tryp_raw_object@meta.data$nFeature_RNA, )) + geom_histogram(binwidth = 20) + scale_x_continuous(name = "Feature count", breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000))  + geom_vline(aes(xintercept=2000), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
ggplot(tryp_raw_object@meta.data, aes(x = tryp_raw_object@meta.data$nCount_RNA, )) + geom_histogram(binwidth = 20) + scale_x_continuous(name = "UMI count", breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000)) + geom_vline(aes(xintercept=3000), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
ggplot(tryp_raw_object@meta.data, aes(x = tryp_raw_object@meta.data$percent.MT, )) + geom_histogram(binwidth = 0.01) + scale_x_log10(name = "% kDNA genes", breaks = c(0.25, 0.5, 1, 2, 3, 5)) + geom_vline(aes(xintercept=2), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
VlnPlot(tryp_raw_object, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)


# Select a subset of cells
tryp <- subset(tryp_raw_object, subset = nFeature_RNA > 50 & nFeature_RNA < 2500 & nCount_RNA < 4000 & percent.MT < 2)

## Export to a SingleCellExperiment object 
sce_tryp <- as.SingleCellExperiment(tryp)

## Normalise data with SCRAN

# Pre-cluster cells. Factors are first generated within clusters then rescaled to normalize between clusters
qclust <- scran::quickCluster(sce_tryp, min.size = 30)

# Compute size factors - removes low abundance genes
sce_tryp <- scran::computeSumFactors(sce_tryp, clusters = qclust)
sce_tryp <- scater::logNormCounts(sce_tryp)

## Convert back to Seurat object 
tryp <- as.Seurat(sce_tryp, counts = "counts", data = "logcounts")

# Find variable features
tryp.features <- FindVariableFeatures(tryp, selection.method = "vst", nfeatures = 2000) 
plot1 <- VariableFeaturePlot(tryp.features, pt.size = 0.1)
top10 <- head(VariableFeatures(tryp.features), 10)
top10_names <- c("glyceraldehyde 3-phosphate dehydrogenase", "flagellum targeting protein kharon1", "cyclophilin a",
                 "Cyb", "pyruvate kinase 1", "ATP-dependent 6-phosphofructokinase",
                 "hypothetical protein, conserved", "GPEET procyclin", "DNA polymerase kappa", "beta tubulin")
top2000 <- head(VariableFeatures(tryp.features), 2000)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10_names, repel = TRUE)

tryp <- FindVariableFeatures(tryp, selection.method = "vst", nfeatures = 2000)

# Scale data and regress the total UMI variation
all.genes <- rownames(tryp)
tryp <- ScaleData(tryp, vars.to.regress = c("nCount_RNA"), do.scale = TRUE, features = all.genes)

# Run PCA using just the top 2000 variable genes
tryp <- RunPCA(tryp, verbose = FALSE, )

# Determine percent of variation associated with each PC
pct <- tryp[["pca"]]@stdev / sum(tryp[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + xlab("Cumulative % variance") + ylab("% of variance") + guides(color=guide_legend(title="> 0.05% variance")) +
  geom_vline(xintercept = 90, color = "black", linetype="dashed") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), 
        axis.text.y =  element_text(size = 10), legend.text = element_text(size = 10))

## Use the appropriate number of dimentions found to find nearest neighbours

tryp <- FindNeighbors(tryp, dims = 1:6, k.param = 30, nn.method = "annoy", annoy.metric = "euclidean")

## Identify the appropriate resolution of clusters
# Generate a copy seurat object and try several resolutions 
tryp_copy <- FindClusters(tryp, resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1))

# Identify which clusters persist at several resolutions with Clustree
clustree(tryp_copy, prefix = "RNA_snn_res.", exprs = c("data", "counts", "scale.data"), assay = NULL)

## Find clusters and produce UMAP
# Pick resolution based on clustree result
tryp <- FindClusters(tryp, resolution = 0.4)
# Run UMAP with same number of dims as for FindNeighbors
tryp <- RunUMAP(tryp, dims = 1:6, reduction = "pca", min.dist = 0.01, n.neighbors = 50, n.epochs = 500, umap.method = "umap-learn")
DimPlot(tryp, reduction = "umap", label = TRUE,
        label.size = 6,
        plot.title = "UMAP") + NoLegend()
DimPlot(tryp, reduction = "umap", label = TRUE,
        label.size = 6,
        plot.title = "UMAP", colo )

FeaturePlot(object = tryp, features = "Tb927.7.5940") #PAD2
FeaturePlot(object = tryp, features = "Tb927.7.2660") #ZC3H20
FeaturePlot(object = tryp, features = "Tb927.6.4280") #GAPDH
VlnPlot(object = tryp, features = "Tb927.6.4280")
FeaturePlot(object = tryp, features = "nCount_RNA")



save(tryp, file = "Seurat_0-6res")
# Find markers
tryp.markers <- FindAllMarkers(tryp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(tryp.markers, file = "res0-6_markers")
top.markers <- tryp.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
levels(tryp) <- c("5", "3", "1", "0", "4", "2")
DoHeatmap(tryp, features = top.markers$gene) + scale_fill_gradientn(colors = c("purple", "white", "yellow"))


# Name clusters
new.cluster.ids <- c("Slender a", "Slender b", "Intermediate", "Stumpy a", "Stumpy b", "low UMI")
names(new.cluster.ids) <- levels(tryp)
tryp <- RenameIdents(tryp, new.cluster.ids)
DimPlot(tryp, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Subset

tryp_subset <- subset(tryp, idents = c("Slender a", "Slender b"))
FeaturePlot(object = tryp_subset, features = "Tb927.8.2780") #PAD2
FeaturePlot(object = tryp_subset, features = "Tb927.7.2660") #ZC3H20
FeaturePlot(object = tryp_subset, features = "Tb927.6.4280") #GAPDH
slenders <- subset(tryp, subset = Tb927.6.4280 > 1.5)
