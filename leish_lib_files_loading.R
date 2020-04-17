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
#library(shinydashboard)
library(shiny)
library(DT)
library(shinythemes)
library(shinycssloaders)



# adapt to your path
setwd("D:/GCRF_UoG/Felix_raw_data/")

##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
leishmania_umap_list_res = readRDS("D:/GCRF_UoG/Felix_raw_data/leishmania_umap_list_res.rds")

##Reading in the list of precomputed table of cluster markers
leishmania_marker_tables_res = readRDS("D:/GCRF_UoG/Felix_raw_data/leishmania_marker_tables_res.rds")

##Reading in all the genes
all_genes = readRDS("D:/GCRF_UoG/Felix_raw_data/all_genes.rds")

##Reading in and processing table of uniprot links for all mouse genes
# uniprot_info_raw = fread("/GCRF_UoG/Felix_raw_data/uniprot table/unipro-mouseID")
# uniprot_info_raw$uniprot = paste('<a href="https://www.uniprot.org/uniprot/',uniprot_info_raw$Entry,'" target="_blank">', uniprot_info_raw$Entry,'</a>', sep = "")
# write.table(uniprot_info_raw, "/GCRF_UoG/Felix_raw_data/uniprot_info_with_link", row.names = F, sep = "\t", quote = F)
# uniprot_info = fread("D:/GCRF_UoG/Felix_raw_data/uniprot_info_with_link", stringsAsFactors = F)

##Declaring and assigning variables
dim=8

res1 = 0.2
res2 = 0.7
diff_res = 0.1
cluster_names = c("unknown","unknown")

# cluster_names = c("Il17a +ve cells", "Ccr7 +ve cells", "Ly6c2 +ve cells", "Gzma +ve cells", "Cdk6 +ve cells")
# fav_genes = c("Cdk6", "Gzma", "Ly6c2", "Ccr7", "Il17a")
# conditions = c("CD18 KO", "WT")

##Code block for processing raw 10X data and computing the Seurat objects and cluster marker tables
#if (!(all(exists("tcells.combined.umap.list"), exists("tcells.combined.clusters.tables")))) {
if (!(all(exists("leishmania_umap_list_res"), exists("leishmania_marker_tables_res"), exists("all_genes")))) {
  
  # Load data set
  leish_raw = Read10X(data.dir = "D:/GCRF_UoG/Felix_raw_data/leish_dset/")
  
  # Initialize the Seurat object with the raw (non-normalized data).
  leish_raw_object = CreateSeuratObject(counts = leish_raw, project = "Cells", min.cells = 5, min.features = 3)
  
  ##Cell death genes
  #mt_marker_genes = c("ND7", "COIII", "MURF4", "MURF1", "ND1", "MURF2", "COI", "COII", "ND4", "RPS12", "ND5", "Cyb")
  mt_marker_genes = c("ND7", "ND4", "CYb")
  leish_raw_object[["percent.MT"]] = PercentageFeatureSet(leish_raw_object, features = mt_marker_genes)
  
  # Visualize QC metrics as a violin plot, need way of using mito genes
  # featureVln = VlnPlot(leish_raw_object, features = "nFeature_RNA", pt.size = 0.01) + labs(title = "Features") + ylab("Count") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # UMIVln = VlnPlot(leish_raw_object, features = "nCount_RNA", pt.size = 0.01) + labs(title = "UMIs") + ylab("Count") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # MTVln = VlnPlot(leish_raw_object, features = "percent.MT", pt.size = 0.01) + labs(title = "kDNA Genes") + ylab("% Features") + xlab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_blank(), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 10), axis.text.y =  element_text(size = 10))
  # ggarrange(featureVln, UMIVln, MTVln,ncol = 3, legend = "none", label.x = 0)
  # 
  # UMI_MT = FeatureScatter(leish_raw_object, feature1 = "nCount_RNA", feature2 = "percent.MT", pt.size = 0.1, smooth = FALSE) + labs(title = "") + scale_y_continuous(name = "% kDNA genes") + xlab("UMI Count") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # UMI_feature = FeatureScatter(leish_raw_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)+ labs(title = "") + ylab("Feature Count") + xlab("UMI Count") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # ggarrange(UMI_MT, UMI_feature,  ncol = 2, legend = "none")
  # 
  # ggplot(leish_raw_object@meta.data, aes(x = leish_raw_object@meta.data$nFeature_RNA, )) + geom_histogram(binwidth = 20) + scale_x_continuous(name = "Feature count", breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000))  + geom_vline(aes(xintercept=2000), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # ggplot(leish_raw_object@meta.data, aes(x = leish_raw_object@meta.data$nCount_RNA, )) + geom_histogram(binwidth = 20) + scale_x_continuous(name = "UMI count", breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000)) + geom_vline(aes(xintercept=3000), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # ggplot(leish_raw_object@meta.data, aes(x = leish_raw_object@meta.data$percent.MT, )) + geom_histogram(binwidth = 0.01) + scale_x_log10(name = "% kDNA genes", breaks = c(0.25, 0.5, 1, 2, 3, 5)) + geom_vline(aes(xintercept=.04), color = "red") + ylab("Cells") + theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), axis.text.y =  element_text(size = 10))
  # VlnPlot(leish_raw_object, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  
  
  # Select a subset of cells
  leish = subset(leish_raw_object, subset = nFeature_RNA > 50 & nFeature_RNA < 2500 & nCount_RNA < 3500 & percent.MT < .04)
  
  ## Export to a SingleCellExperiment object 
  sce_leish = as.SingleCellExperiment(leish)
  
  ## Normalise data with SCRAN
  
  # Pre-cluster cells. Factors are first generated within clusters then rescaled to normalize between clusters
  qclust = scran::quickCluster(sce_leish, min.size = 30)
  
  # Compute size factors - removes low abundance genes
  sce_leish = scran::computeSumFactors(sce_leish, clusters = qclust)
  sce_leish = scater::logNormCounts(sce_leish)
  
  ## Convert back to Seurat object 
  leish = as.Seurat(sce_leish, counts = "counts", data = "logcounts")
  
  # Find variable features
  leish.features = FindVariableFeatures(leish, selection.method = "vst", nfeatures = 2000) 
  # plot1 = VariableFeaturePlot(leish.features, pt.size = 0.1)
  # top10 = head(VariableFeatures(leish.features), 10)
  # 
  # # top10_names = c("glyceraldehyde 3-phosphate dehydrogenase", "flagellum targeting protein kharon1", "cyclophilin a",
  # #                  "Cyb", "pyruvate kinase 1", "ATP-dependent 6-phosphofructokinase",
  # #                  "hypothetical protein, conserved", "GPEET procyclin", "DNA polymerase kappa", "beta tubulin")
  # 
  # # top2000 = head(VariableFeatures(leish.features), 2000)
  # plot1
  # # plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  leish = FindVariableFeatures(leish, selection.method = "vst", nfeatures = 2000)
  # Scale data and regress the total UMI variation
  
  ##capturing all genes at this stage
  all_genes = rownames(leish)
  saveRDS(all_genes, "all_genes.rds")
  
  leish = ScaleData(leish, vars.to.regress = c("nCount_RNA"), do.scale = TRUE, features = all_genes)
  
  # Run PCA using just the top 2000 variable genes
  leish = RunPCA(leish, verbose = FALSE)
  
  # # Determine percent of variation associated with each PC
  # pct = leish[["pca"]]@stdev / sum(leish[["pca"]]@stdev) * 100
  # # Calculate cumulative percents for each PC
  # cumu = cumsum(pct)
  # # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  # co1 = which(cumu > 90 & pct < 5)[1]
  # co1
  # # Determine the difference between variation of PC and subsequent PC
  # co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
  # # last point where change of % of variation is more than 0.1%.
  # co2
  # # Minimum of the two calculation
  # pcs = min(co1, co2)
  # pcs
  # # Create a dataframe with values
  # plot_df = data.frame(pct = pct, 
  #                       cumu = cumu, 
  #                       rank = 1:length(pct))
  # # Elbow plot to visualize 
  # ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  #   geom_text() + xlab("Cumulative % variance") + ylab("% of variance") + guides(color=guide_legend(title="> 0.05% variance")) +
  #   geom_vline(xintercept = 90, color = "black", linetype="dashed") + 
  #   geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  #   theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10),
  #         axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), 
  #         axis.text.y =  element_text(size = 10), legend.text = element_text(size = 10))
  # 
  # ## Use the appropriate number of dimentions found to find nearest neighbours
  
  leish = FindNeighbors(leish, dims = 1:dim, k.param = 30, nn.method = "annoy", annoy.metric = "euclidean")
  
  ## Identify the appropriate resolution of clusters
  # Generate a copy seurat object and try several resolutions 
  #leish_copy = FindClusters(leish, resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1))
  
  # Identify which clusters persist at several resolutions with Clustree
  #clustree(leish_copy, prefix = "RNA_snn_res.", exprs = c("data", "counts", "scale.data"), assay = NULL)
  
  ## Find clusters and produce UMAP
  # Pick resolution based on clustree result
  #leish = FindClusters(leish, resolution = 0.6)
  
  ##Precomputing and saving the list of Seurat objects with different clusters through adjusting of resolution from 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 & 0.8
  leishmania_umap_list_res = lapply(seq(res1, res2, by = diff_res), function(x) {
    leish_res = FindClusters(leish, resolution = x)
    RunUMAP(leish_res, dims = 1:dim, reduction = "pca", min.dist = 0.01, n.neighbors = 50, n.epochs = 500, umap.method = "umap-learn")
    
  })
  saveRDS(leishmania_umap_list_res, "D:/GCRF_UoG/Felix_raw_data/leishmania_umap_list_res.rds")
  
  
  ##Precomputing and saving the list of tables of cluster markers
  leishmania_marker_tables_res = lapply(leishmania_umap_list_res, function(x) { 
    #lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      leish.markers = FindAllMarkers(x, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
      
    #})
  })
  
  saveRDS(leishmania_marker_tables_res, "D:/GCRF_UoG/Felix_raw_data/leishmania_marker_tables_res.rds")
  
  
  # Run UMAP with same number of dims as for FindNeighbors
  # leish = RunUMAP(leish, dims = 1:6, reduction = "pca", min.dist = 0.01, n.neighbors = 50, n.epochs = 500, umap.method = "umap-learn")
  # DimPlot(leish, reduction = "umap", label = TRUE,
  #         label.size = 6,
  #         plot.title = "UMAP") + NoLegend()
  # DimPlot(leish, reduction = "umap", label = TRUE,
  #         label.size = 6,
  #         plot.title = "UMAP")
  # FeaturePlot(object = leish, features = "nCount_RNA")
  # 
  # leish.markers = FindAllMarkers(leish, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
  # write.csv(leish.markers, file = "res0-6_markers")
  # top.markers = leish.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # top10 = leish.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # levels(leish) = c("5", "3", "1", "0", "4", "2")
  # DoHeatmap(leish, features = top.markers$gene) + scale_fill_gradientn(colors = c("purple", "white", "yellow"))
  # 
  
  # Name clusters
  # names(new.cluster.ids) = levels(leish)
  # leish = RenameIdents(leish, new.cluster.ids)
  # DimPlot(leish, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  
  
}

which_numeric_cols = function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}