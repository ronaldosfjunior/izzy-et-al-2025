# calling packages
library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(Matrix)
library(celldex)
library(scater)
library(ensembldb)
library(car)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(12))
library(harmony)
#library(SeuratDisk)

setwd("/media/patrick/GERVAZIO/Bioinfo/weiner_lab/izzy/scRNAseq/mice/Seurat/")

file.dir <- "../h5_files/"
files.set <- c(
  "ADBl1Tg1",
  "ADBl1Tg2",
  "ADBl1Tg3",
  "ADBl1Tg4",
  "ADBl1wt1",
  "ADBl1wt2",
  "ADBl1wt3",
  "ADBl1wt4"
)

data.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i], "_raw_feature_bc_matrix.h5")
  data <- Read10X_h5(filename = path)
  dataset_name <- files.set[i]
  condition <- sub("ADBl1(.*)[0-9]", "\\1", files.set[i])
  group <- sub("ADBl1.*([0-9])", "\\1", files.set[i])
  mat <- CreateSeuratObject(counts = data, min.cells = 5, min.features = 500) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 20) %>%
    #NormalizeData(normalization.method = "LogNormalize") %>%
    #FindVariableFeatures(selection.method = "vst")
    SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
  #RunPCA(npcs = 50)
  mat$condition <- condition
  mat$patient <- dataset_name
  mat$group <- group
  data.list[[i]] <- mat
  rm(mat)
  rm(data)
}
# Name layers
names(data.list) <- files.set

# Merge layers and create seurat obj during merging 
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 2000)
prot.combined <- merge(data.list[[1]], y = data.list[2:length(data.list)], 
                       add.cell.ids = files.set, merge.data = T)

VariableFeatures(prot.combined) <- features

# Normalize and scale merged obj
prot.combined <- NormalizeData(prot.combined)
#prot.combined <- FindVariableFeatures(prot.combined)
#prot.combined <- ScaleData(prot.combined)

## Dimensionality reduction and integration
prot.combined <- RunPCA(prot.combined)
gc()
prot.combined <- FindNeighbors(prot.combined, dims = 1:30, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 0.5, cluster.name = "unintegrated_clusters")
prot.combined <- RunUMAP(prot.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()
DimPlot(prot.combined, reduction = "umap.unintegrated", label = T, raster=F,
        #group.by = "condition"
        split.by = "patient", ncol = 4
)
gc()

prot.combined <- RunHarmony(prot.combined, group.by.vars = "patient")
#prot.combined <- RunHarmony(prot.combined, group.by.vars = "patient", #assay.use = "SCT", reduction.save = "harmony.umap", reduction = "umap.unintegrated")
#prot.combined <- RunHarmony(prot.combined, group.by.vars = "patient", assay.use = "SCT", reduction.save = "harmony.SCT")
#prot.combined <- RunHarmony(prot.combined, group.by.vars = "patient", assay.use = "SCT", reduction.save = "harmony.SCT.umap", reduction = "umap.unintegrated")
prot.combined <- RunUMAP(prot.combined, dims = 1:30, reduction = "harmony", reduction.name = "umap")
prot.combined <- FindNeighbors(prot.combined, dims = 1:2, reduction = "umap")
prot.combined <- FindClusters(prot.combined, resolution = 0.3, cluster.name = "harmony_clusters")
gc()
DimPlot(prot.combined, reduction = "umap", label = T, raster=F,
        #group.by = "condition"
        split.by = "condition"#, ncol = 4
)

############## Identify cell types with scType

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = prot.combined[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results = do.call("rbind", lapply(unique(prot.combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(prot.combined@meta.data[prot.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(prot.combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

prot.combined@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  prot.combined@meta.data$customclassif[prot.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Visualization on UMAP
DimPlot(prot.combined, reduction = "umap", label = TRUE, repel = TRUE, 
        group.by = 'customclassif', 
        #split.by = "patient"
)   
prot.combined2 <- RenameIdents(prot.combined, `0` = "Classical Mono",
                               `2` = "Classical Mono",
                               `20` = "Intermediate Mono",
                               `18` = "Nonclassical Mono"
)


## save the Seurat object
saveRDS(prot.combined, file = "./AD-mice.izzy.seurat-harmony.rds")
#DefaultAssay(prot.combined) <- "RNA"
#SaveH5Seurat(prot.combined, filename = "AD.izzy.seurat-harmony.h5Seurat")
#Convert("AD.izzy.seurat-harmony.h5Seurat", dest = "h5ad")
## load back Seurat object
prot.combined <- readRDS("AD-mice.izzy.seurat-harmony.rds")

prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = T, min.pct = 0.25, logfc.threshold = 0.33) %>% group_by(cluster)
write.csv(as.data.frame(prot.markers), file="All_pos_markers.csv")


VlnPlot(prot.combined, features = c("Cd14","Fcgr3"), 
        #split.by = "condition",
        #idents = c("Classical Mono","Intermediate Mono","Nonclassical Mono")
        )#assay = "SCT", 
FeaturePlot(prot.combined, features = c("Cd14","Fcgr3"), reduction = "umap")#, split.by = "condition")

DimPlot(prot.combined, reduction = "umap", label = T, raster=F,
        #group.by = "condition"
        split.by = "patient", ncol = 4
)

prot.response <- FindMarkers(prot.combined, ident.1 = "20_AD", 
                             ident.2 = "20_C",
                             slot = "data",
                             assay = "RNA",
                             features = NULL,
                             logfc.threshold = 0,
                             test.use = "wilcox",
                             min.pct = 0.0,
                             min.diff.pct = -Inf,
                             verbose = TRUE,
                             only.pos = FALSE,
                             max.cells.per.ident = Inf,
                             random.seed = 1,
                             latent.vars = NULL,
                             min.cells.feature = 3,
                             min.cells.group = 3,
                             pseudocount.use = 1,
                             mean.fxn = NULL,
                             fc.name = NULL,
                             base = 2,
                             densify = FALSE,
                             recorrect_umi = TRUE)

write.csv(as.data.frame(prot.response), file="ADxC_intermediate.mono_DEGs.csv")


prot.response <- FindMarkers(prot.combined, ident.1 = "18_AD", 
                             ident.2 = "18_C",
                             slot = "data",
                             assay = "RNA",
                             features = NULL,
                             logfc.threshold = 0,
                             test.use = "wilcox",
                             min.pct = 0.0,
                             min.diff.pct = -Inf,
                             verbose = TRUE,
                             only.pos = FALSE,
                             max.cells.per.ident = Inf,
                             random.seed = 1,
                             latent.vars = NULL,
                             min.cells.feature = 3,
                             min.cells.group = 3,
                             pseudocount.use = 1,
                             mean.fxn = NULL,
                             fc.name = NULL,
                             base = 2,
                             densify = FALSE,
                             recorrect_umi = TRUE)

write.csv(as.data.frame(prot.response), file="ADxC_nonclassical.mono_DEGs.csv")

prot.response <- FindMarkers(prot.combined, ident.1 = c("0_AD","2_AD"), 
                             ident.2 = c("0_C","2_C"),
                             slot = "data",
                             assay = "RNA",
                             features = NULL,
                             logfc.threshold = 0,
                             test.use = "wilcox",
                             min.pct = 0.0,
                             min.diff.pct = -Inf,
                             verbose = TRUE,
                             only.pos = FALSE,
                             max.cells.per.ident = Inf,
                             random.seed = 1,
                             latent.vars = NULL,
                             min.cells.feature = 3,
                             min.cells.group = 3,
                             pseudocount.use = 1,
                             mean.fxn = NULL,
                             fc.name = NULL,
                             base = 2,
                             densify = FALSE,
                             recorrect_umi = TRUE)

write.csv(as.data.frame(prot.response), file="ADxC_classical.mono_DEGs.csv")








zk.combined2$celltypes <- sub("(.*)_.*", "\\1", Idents(zk.combined2))
zk.combined2$celltypes.condition <- paste(zk.combined2$celltypes, zk.combined2$condition, sep = "_")
zk.combined2$celltypes.condition <- paste(zk.combined2$customclassif, zk.combined2$condition, sep = "_")
prot.combined$cluster.condition <- paste(Idents(prot.combined), prot.combined$condition, sep = "_")
zk.combined2$cluster.condition <- paste(zk.combined2$seurat_clusters, zk.combined2$condition, sep = "_")
Idents(prot.combined) <- "seurat_clusters"
Idents(zk.combined2) <- "celltypes"
Idents(zk.combined2) <- "celltypes.condition"
Idents(prot.combined) <- "cluster.condition"
Idents(zk.combined2) <- "customclassif"




