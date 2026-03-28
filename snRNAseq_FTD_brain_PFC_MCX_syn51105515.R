# Seurat v5
library(openxlsx)
library(MAST)
library(Seurat)
library(harmony)
library(reticulate)
library(BPCells)
library(dplyr)
library(RColorBrewer)
library(SeuratWrappers)
library(Azimuth)
library(Matrix)
library(sctransform)
library(car)
library(scater)
library(ggplot2)
library(patchwork)
library(ggrepel)
options(future.globals.maxSize = 3e+09)
options(Seurat.object.assay.version = "v5")
library(BiocParallel)
register(MulticoreParam(14))
library(paletteer)
library(nichenetr)
library(tidyverse)
library(circlize)
library(ggpubr)
library(VennDiagram)
library(pheatmap)
library(gplots)
library(GOplot)
library(readr)
library(robustbase)
library(EnhancedVolcano)
library(alluvial)
library(reshape2)
library(scico)
library(ggsci)
library(rcartocolor)
library(ggside)
library(viridis)
library(ggstatsplot)
library(grid)
library(shadowtext)
library(tidyr)


setwd("/media/patrick/JANDERSON/Bioinfo/weiner_lab/public_data/ALS/brain/Seurat/")


file.dir <- "../processed/"
data.list <- c()
meta1 <- read.delim("meta_FTD.txt", header = T)

files.set <- meta1$Sample_ID

# files.set <- t(read.delim("FTLD_PN.txt", header = F))

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i], "/")
  data <- Read10X(data.dir = path)
  dataset_name <- files.set[i]
  # condition <- sub("(.*)-.*", "\\1", files.set[i])
  # treatment1 <- treatment[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 500) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5) %>%
    NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst")
  #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
  #RunPCA(npcs = 50)
  # mat$treatment <- treatment1
  mat$sample <- dataset_name
  data.list[[i]] <- mat
  rm(mat)
  rm(data)
}
# Name layers
names(data.list) <- files.set

# Merge layers and create seurat obj during merging 
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
prot.combined <- merge(data.list[[1]], y = data.list[2:length(data.list)], 
                       add.cell.ids = files.set, merge.data = T)


rm(data.list)
gc()
VariableFeatures(prot.combined) <- features

# Visualize QC metrics as a violin plot
#VlnPlot(prot.combined, pt.size = 0.0,group.by = "cohort",raster = F,
#        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",raster = F)
#plot2 <- FeatureScatter(prot.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster = F)
#plot1 + plot2

### Normalize and scale merged obj
prot.combined <- NormalizeData(prot.combined, normalization.method = "LogNormalize")
prot.combined <- FindVariableFeatures(prot.combined, selection.method = "vst", nfeatures = 3000)
prot.combined <- ScaleData(prot.combined, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                           model.use = "linear")#, features = all.genes)
gc()
### Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(prot.combined), 10)

### plot variable features with and without labels
#plot1 <- VariableFeaturePlot(prot.combined, raster = F)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, raster = F)
#plot1 + plot2

### Dimensionality reduction and integration
prot.combined <- RunPCA(prot.combined, npcs = 50)
gc()
ElbowPlot(prot.combined, ndims = 50)
prot.combined <- FindNeighbors(prot.combined, dims = 1:40, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 2.5, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        ncol = 2,
        label = T,
        repel = T,
        #group.by = "seurat_clusters"
        # split.by = "disease"
)#, combine = F)

### add metadata
prot.combined$disease <- prot.combined$sample 
prot.combined$sex <- prot.combined$sample
prot.combined$disease2 <- prot.combined$sample
prot.combined$donor <- prot.combined$sample
prot.combined$region <- prot.combined$sample

for (i in 1:length(meta1$Sample_ID)) {
  prot.combined$disease <- recode(prot.combined$disease, "meta1$Sample_ID[i] = meta1$condition[i]")
  prot.combined$sex <- recode(prot.combined$sex, "meta1$Sample_ID[i] = meta1$sex[i]")
  prot.combined$disease2 <- recode(prot.combined$disease2, "meta1$Sample_ID[i] = meta1$group[i]")
  prot.combined$donor <- recode(prot.combined$donor, "meta1$Sample_ID[i] = meta1$Donor[i]")
  prot.combined$region <- recode(prot.combined$region, "meta1$Sample_ID[i] = meta1$Region[i]")
}
prot.combined$region_disease <- paste(prot.combined$region, prot.combined$disease, sep = "_")


## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats_FTD_PFC_MCX.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats_FTD_PFC_MCX.Rds")


## cell type annotation
prot.combined[["RNA"]] <- JoinLayers(prot.combined[["RNA"]])
prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap", reorder = T)
prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = T, min.pct = 0.35, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.clusters.markers.xlsx")
prot.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5


# DotPlot(prot.combined, features = unique(top5$gene),
#         col.max = 20, 
#         dot.scale = 10, 
#         #cluster.idents = F, #group.by = "patho",
#         #scale = F,
#         #split.by = "disease"
# ) + RotatedAxis() + coord_flip()

cell.types <- c(
  "CD8+ T cells",
  "Monocytes",
  "DAM Microglia",
  "Homeostatic Microglia",
  "Astrocytes",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Vascular smooth muscle cells",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Fibroblasts"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = cell.types)

DotPlot(prot.combined, features = c( 
                                     "FBLN1","FBLN5", # fibroblasts
                                     "GAD1","GAD2",#"TAC1","PENK","SST","NPY","MYBPC1","PVALB","GABBR2",  "SLC32A1",  #Inhibitory_neurons
                                     "SLC17A7","NRGN","CAMK2A","HTR2C", #"SATB2", #"COL5A1","SDK2","NEFM", "SLC17A6",   #Excitatory_neurons
                                     "TAGLN","MYH11", # vascular smooth muscle cells
                                     "PDE5A","PTH1R","P2RY14",#"AMBP","HIGD1B","COX4I2", "AOC3","ABCC9","KCNJ8","CD248",  #Pericytes
                                     "FLT1","CLDN5", #"VTN","ITM2A", "VWF", "FAM167B","BMX","CLEC1B",    #Endothelial_cells
                                     "VCAN","PCDH15",#"CSPG4","PDGFRA", "SOX10","NEU4", "GPR37L1","C1QL1","CDO1","EPN2",   #Oligodendrocyte_precursor_cells
                                     "MBP","MOBP","PLP1",#"MOG","CLDN11","MYRF","GALC","ERMN","MAG", "OLIG2",   #Oligodendrocytes
                                     "GFAP","AQP4",#"LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4", "SLC1A3",  #Astrocytes
                                     "P2RY12","CSF1R","C3","APOE", #,"CD74","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5","AIF1","IL1B","IRF8",   #Microglia
                                     "MS4A4A","CD163","MAFB","TNFAIP2","IL15","ASAH1","PLA2G7", #"PLXND1","EMILIN2","SIGLEC1","F13A1","MARCO","GAS7", # monocytes
                                     #"CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2"#, #Macrophages
                                     "PTPRC","PRF1","CD8A","CD3E","CD247" # T cells
                                     #"MMP9", # neutrophil
                                     # "CD8A", "CD3E",
                                     
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = F, group.by = "celltypes",
#scale = F,
#split.by = "cohort"
) + RotatedAxis() + coord_flip()

DotPlot(prot.combined, features = c("P2RY12","CSF1R","CD74","C3",#"CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",
                               "IL1B","APOE","IRF8", #"CST7","CLEC7A","AIF1",
                               "MARCO","MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2","CD4", #"GDA",
                               "CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7"#,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
),
col.max = 20,
dot.scale = 10, 
cluster.idents = T, #group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis() #+ coord_flip()

Idents(prot.combined) <- "seurat_clusters"
prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Astrocytes_1",
                              `1`= "Oligodendrocytes_1",
                              `2`= "Oligodendrocytes_2",
                              `3`= "Oligodendrocytes precursor cells_1",
                              `4`= "Oligodendrocytes_3",
                              `5`= "Inhibitory neurons_1",
                              `6`= "Oligodendrocytes_4",
                              `7`= "Excitatory neurons_1",
                              `8`= "Excitatory neurons_2",
                              `9`= "Oligodendrocytes precursor cells_2",
                              `10`= "Oligodendrocytes_5",
                              `11`= "Oligodendrocytes_6",
                              `12`= "Excitatory neurons_3",
                              `13`= "Astrocytes_2",
                              `14`= "Astrocytes_3",
                              `15`= "Excitatory neurons_4",
                              `16`= "Inhibitory neurons_2",
                              `17`= "Homeostatic Microglia_1",
                              `18`= "Inhibitory neurons_3",
                              `19`= "Oligodendrocytes precursor cells_3",
                              `20`= "Excitatory neurons_5",
                              `21`= "DAM Microglia_1",
                              `22`= "Inhibitory neurons_4",
                              `23`= "Inhibitory neurons_5",
                              `24`= "Excitatory neurons_6",
                              `25`= "Oligodendrocytes_7",
                              `26`= "Excitatory neurons_7",
                              `27`= "Inhibitory neurons_6",
                              `28`= "Inhibitory neurons_7",
                              `29`= "Excitatory neurons_8",
                              `30`= "Astrocytes_4",
                              `31`= "Homeostatic Microglia_2",
                              `32`= "Monocytes_1",
                              `33`= "Inhibitory neurons_8",          ### 
                              `34`= "Endothelial cells_1",
                              `35`= "Astrocytes_5",
                              `36`= "Oligodendrocytes precursor cells_4",
                              `37`= "Astrocytes_6",
                              `38`= "Astrocytes_7",              #### 
                              `39`= "Inhibitory neurons_9",
                              `40`= "Excitatory neurons_9",
                              `41`= "Excitatory neurons_10",       #########
                              `42`= "Oligodendrocytes precursor cells_5",
                              `43`= "Inhibitory neurons_10",
                              `44`= "Excitatory neurons_11",       
                              `45`= "Inhibitory neurons_11",
                              `46`= "Homeostatic Microglia_3",
                              `47`= "Astrocytes_8",
                              `48`= "Excitatory neurons_12",
                              `49`= "Pericytes_1",          ######
                              `50`= "Excitatory neurons_13",
                              `51`= "Astrocytes_9",
                              `52`= "Inhibitory neurons_12",
                              `53`= "Fibroblasts_1",
                              `54`= "Excitatory neurons_14",
                              `55`= "Excitatory neurons_15",
                              `56`= "Inhibitory neurons_13",
                              `57`= "Neuroendocrine cells_1",
                              `58`= "Excitatory neurons_16",
                              `59`= "Vascular smooth muscle cells_1",
                              `60`= "Oligodendrocytes precursor cells_6",
                              `61`= "Astrocytes_10",
                              `62`= "CD8+ T cells_1"
                              )

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "seurat_clusters"

## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats_FTD_PFC_MCX_annotated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats_FTD_PFC_MCX_annotated.Rds")

################## IL15 expression by cell types
# Idents(prot.combined) <- "cohort"
gene1 <- "IL15"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "celltypes", #cols = colors,
              # idents = "Hong Kong University"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.000001, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1


FeaturePlot(prot.combined, features = gene1, raster = T, pt.size = 1.6,
            cols = c("gray80","#9E1021"),
            # split.by = "cohort",
            # ncol = 2,
            # min.cutoff = 0.5,
            # max.cutoff = 2.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")


#####################################


############# Fig 6 

celltypes <- c(
  "CD8+ T cells",
  "Monocytes",
  "DAM Microglia",
  "Homeostatic Microglia",
  "Astrocytes",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Vascular smooth muscle cells",
  "Fibroblasts",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = celltypes)

col1 <- c(
  paletteer_c("ggthemes::Red-Green-Gold Diverging", 14)[14:1] 
)

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = T, pt.size = 0.9, 
        cols = col1,
        # ncol = 2,
        label = T,
        repel = T,
        group.by = "celltypes",
        # split.by = "region_disease"
)

##### DEGs analysis

cell.types <- c(
  "CD8+ T cells",
  "Monocytes",
  "DAM Microglia",
  "Homeostatic Microglia",
  "Astrocytes",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Vascular smooth muscle cells",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons",
  "Fibroblasts"
)

prot.combined$celltypes.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "celltypes.disease"

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_FTLD"), 
                              ident.2 = paste0(cell.types[i], "_PN"),
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
                              recorrect_umi = TRUE
  )
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_MCX_PFC_FTLDxPN_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}


Idents(prot.combined) <- "seurat_clusters"
cell.types <- as.character(c(0:62))

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "cluster.disease"

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_FTLD"), 
                              ident.2 = paste0(cell.types[i], "_PN"),
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
                              recorrect_umi = TRUE
  )
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_MCX_PFC_FTLDxPN_clus", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

## PFC only
cell.types <- c(
  "Astrocytes",
  "Homeostatic Microglia",
  "Monocytes",
  "DAM Microglia",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Vascular smooth muscle cells",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons",
  "CD8+ T cells",
  "Fibroblasts"
)

prot.combined$region.celltypes.disease <- paste(prot.combined$celltypes.disease, prot.combined$region, sep = "_")
Idents(prot.combined) <- "region.celltypes.disease"

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_FTLD_MCX"), 
                              ident.2 = paste0(cell.types[i], "_PN_MCX"),
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
                              recorrect_umi = TRUE
  )
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_MCX_FTLDxPN_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

zk.response0 <- FindMarkers(prot.combined, ident.1 = c("Monocytes_FTLD_PFC","DAM Microglia_FTLD_PFC","Homeostatic Microglia_FTLD_PFC"), 
                            ident.2 = c("Monocytes_PN_PFC","DAM Microglia_PN_PFC","Homeostatic Microglia_PN_PFC"),
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
                            recorrect_umi = TRUE
)
#zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
#write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_PFC_FTLDxPN_myeloids_DEGs.xlsx")
rm(zk.response0)


Idents(prot.combined) <- "seurat_clusters"
cell.types <- as.character(c(48:62))

prot.combined$region.cluster.disease <- paste(prot.combined$cluster.disease, prot.combined$region, sep = "_")
Idents(prot.combined) <- "region.cluster.disease"

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_FTLD_MCX"), 
                              ident.2 = paste0(cell.types[i], "_PN_MCX"),
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
                              recorrect_umi = TRUE
  )
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_MCX_FTLDxPN_clus", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

########       DEGs

cell.types <- c(
  "Homeostatic Microglia",
  "Monocytes",
  "DAM Microglia",
  "MgND Microglia"
)

myeloids$celltype.region.disease <- paste(myeloids$celltypes, myeloids$region, myeloids$clinical_dx, sep = "_")
Idents(myeloids) <- "celltype.region.disease"


# cell.types <- as.character(c(0:27))
# 
# myeloids$cluster.region.disease <- paste(myeloids$seurat_clusters, myeloids$region, myeloids$clinical_dx, sep = "_")
# Idents(myeloids) <- "cluster.region.disease"

disease <- c("bvFTD","AD","PSP-S")
region <- c("preCG","insula")

for (i in 1:length(cell.types)) {
  for (e in 1:length(region)) {
    for (x in 1:length(disease)) {
      zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_",region[e],"_",disease[x]), 
                                  ident.2 = paste0(cell.types[i],"_",region[e],"_Control"),
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
                                  recorrect_umi = TRUE
      )
      #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
      #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
      write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_",region[e],"_",disease[x],"xC_clus", cell.types[i], "_DEGs.xlsx"))
      rm(zk.response0)
    }
  }
}

myeloids$celltypes.disease <- paste(myeloids$celltypes, myeloids$clinical_dx, sep = "_")
Idents(myeloids) <- "celltypes.disease"

# myeloids$cluster.disease <- paste(myeloids$seurat_clusters, myeloids$clinical_dx, sep = "_")
# Idents(myeloids) <- "cluster.disease"

for (i in 1:length(cell.types)) {
  for (x in 1:length(disease)) {
    zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(cell.types[i], "_",disease[x]), 
                                ident.2 = paste0(cell.types[i],"_Control"),
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
                                recorrect_umi = TRUE
    )
    #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_myeloids_",disease[x],"xC_clus_", cell.types[i], "_DEGs.xlsx"))
    rm(zk.response0)
  }
}








DotPlot(prot.combined, features = c("KLRB1","TNF","IFNG",
 #"GZMH","GSDMB","IL2RB","LTB","PDCD1","FASLG","BATF", "IFNG-AS1", "IL2RG","PRF1","GNLY","NKG7","KLRK1","CD3E",
  "CD8A","GZMA","GZMB","TYROBP","CD8B", "IL7R","TOX","GZMK","IRF2","RUNX2",
  "NCK2",
  "CD2",
  "RHOA",
  "ITK",
  "NFATC2",
  "IKBKB",
  "CD96",
  "LAMP2",
  "FCMR"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "region_disease",
idents = c("CD8+ T cells"),
scale = T,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

prot.combined$disease <- factor(prot.combined$disease, levels = c("PN","FTLD"))

DotPlot(prot.combined, features = c("CD8B","KLRB1","TNF","IFNG",
                                    "IL7R","TOX","GZMK",#"GZMH","GSDMB","IL2RB","LTB","PDCD1","FASLG","BATF", "IFNG-AS1", "IL2RG","PRF1","GNLY","KLRK1","CD3E","NKG7",
                                    "CD8A","GZMB","TYROBP","IRF2","RUNX2", # "GZMA",
                                    "NCK2","CD2","RHOA","ITK",  "NFATC2", "IKBKB","CD96","LAMP2","FCMR"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "disease",
idents = c("CD8+ T cells"),
scale = T,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

prot.combined$region_disease <- factor(prot.combined$region_disease, levels = c("PFC_PN","PFC_FTLD","MCX_PN","MCX_FTLD"))

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2","IRF3","IRF1",
            "PIAS1","PPARA","PPARG","CD36","TLR2","TLR4","MYD88","TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",#,"CD86"
            "IL15","IL15RA")

genes1 <- c("IFNGR1","JAK1","JAK2","IRF1",  ### selection for DAM microglia by region and disease
            "PIAS1","PPARA","PPARG","CD36","TLR4","MYD88","TLR6",
            "MAPK14","MAP2K1","MAP2K6","MAP3K3",  #"MAP2K4","MAP2K5","IRF2","MAP3K5","MAP4K3",
            "MAP4K5","MAPK8", ## "MAPK1","MAP4K4","MAP3K7","MAP3K14","MAP3K20","MAP3K4","MAP3K1",
            "CREB1","CREBBP",#,"CD86","TLR2","IL15RA","IFNGR2","IRF3","STAT1","STAT2","STAT3","LRP1","SREBF1","SREBF2",
            "IL15")

genes1 <- c("JAK1","IRF3",  ### selection for monocytes by region and disease
            "PPARA","PPARG","CD36","TLR4","MYD88","TLR6","SREBF2",
            "MAPK14","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K5", 
            "MAP3K14","MAP3K20","MAP4K3","MAP4K5",  # "MAPK1","MAPK8","MAP4K4","MAP3K7","MAP3K4","MAP2K1",
            "CREB1","CREBBP",#,"CD86","IL15RA","LRP1","SREBF1","TLR2","PIAS1","IRF1", "IRF2","JAK2","STAT1","STAT2","STAT3","IFNGR1","IFNGR2",
            "IL15")

genes1 <- c("JAK1","IRF2","IRF3",  ### selection for homeostatic microglia by region and disease
            "PPARA","PPARG","TLR4","MYD88","TLR6","SREBF2",
            "MAPK14","MAP2K6","MAP3K3","MAP3K4", 
            "MAP4K3","MAP4K5", # "MAPK1","MAPK8","MAP3K20","MAP4K4","MAP3K5","MAP3K7","MAP2K1",
            "CREB1","CREBBP",#,"CD86","IL15RA","LRP1","SREBF1","TLR2","IRF1","PIAS1","STAT1","STAT2","STAT3","IFNGR2","MAP3K14","MAP3K1","CD36","IFNGR1","JAK2",
            "IL15")
DotPlot(prot.combined, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "region_disease",
        idents = c("Monocytes"),
        scale = T,
        # split.by = "region_disease"
) + RotatedAxis() + coord_flip()


FeaturePlot(prot.combined, raster = F, #pt.size = 1.5,
            features = "PRF1", 
            split.by = "region_disease",
            # min.cutoff = 1,
            # max.cutoff = 2.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")

my_comparisons <- list(c("PFC_FTLD","PFC_PN"),c("MCX_FTLD","MCX_PN"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,3,2,3)]

p1 <- VlnPlot(prot.combined, features = "IL15",
              pt.size = 0.05, raster = F, group.by = "region_disease", cols = colors,
              idents = c("Monocytes")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IL15"]])+0.))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1



Idents(prot.combined) <- "orig.ident"
mono.markers <- c("MS4A4A","CD163","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2"#,"CD4","GDA","MARCO","SIGLEC1"
)
DAM1 <- t(read.delim("DAM1.txt", header = F))
DAM2 <- t(read.delim("DAM2.txt", header = F))
MGND <- t(read.delim("MGND.txt", header = F))

prot.combined <- AddModuleScore(
  prot.combined,
  features = mono.markers,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Mono.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

prot.combined <- AddModuleScore(
  prot.combined,
  features = unique(DAM1),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "dam1.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

prot.combined <- AddModuleScore(
  prot.combined,
  features = unique(DAM2),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "dam2.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

prot.combined <- AddModuleScore(
  prot.combined,
  features = unique(MGND),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "MGND.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

FeaturePlot(prot.combined, features = "Mono.markers1", raster = T, pt.size = 1.,
            cols = c("gray80","#9E1021"),
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 1,
            max.cutoff = 2.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")


monocytes <- subset(prot.combined2, subset = Mono.markers1 > 0)

DotPlot(monocytes, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "clinical_dx",
        #idents = c("microglia"),
        scale = T,
        # split.by = "clinical_dx"
) + RotatedAxis() + coord_flip()


################## figures related to AD

prot.combined$celltype.disease.region <- paste(prot.combined$celltypes, prot.combined$disease, prot.combined$region, sep = "_")
Idents(prot.combined) <- "celltype.disease.region"
prot.combined$celltype.disease.region <- factor(prot.combined$celltype.disease.region, levels = c("Homeostatic Microglia_PN_PFC","Homeostatic Microglia_FTLD_PFC",
                                                                                                  "DAM Microglia_PN_PFC","DAM Microglia_FTLD_PFC",
                                                                                                  "Monocytes_PN_PFC","Monocytes_FTLD_PFC",
                                                                                                  "Homeostatic Microglia_PN_MCX","Homeostatic Microglia_FTLD_MCX",
                                                                                                  "DAM Microglia_PN_MCX","DAM Microglia_FTLD_MCX",
                                                                                                  "Monocytes_PN_MCX","Monocytes_FTLD_MCX"))
prot.combined$celltype.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "celltype.disease"
prot.combined$celltype.disease <- factor(prot.combined$celltype.disease, levels = c("Homeostatic Microglia_PN","Homeostatic Microglia_FTLD",
                                                                                                  "DAM Microglia_PN","DAM Microglia_FTLD",
                                                                                                  "Monocytes_PN","Monocytes_FTLD"))

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2","IRF3","IRF1",
            "PIAS1","PPARA","PPARG","CD36","TLR2","TLR4","MYD88","TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",#,"CD86"
            "IL15","IL15RA")

genes1 <- c("JAK1","IRF3",  ### selection for monocytes by region and disease
            "PPARA","PPARG","CD36","TLR4","MYD88","TLR6","SREBF2",
            "MAPK14","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K5", 
            "MAP3K14","MAP3K20","MAP4K3","MAP4K5",  # "MAPK1","MAPK8","MAP4K4","MAP3K7","MAP3K4","MAP2K1",
            "CREB1","CREBBP",#,"CD86","IL15RA","LRP1","SREBF1","TLR2","PIAS1","IRF1", "IRF2","JAK2","STAT1","STAT2","STAT3","IFNGR1","IFNGR2",
            "IL15")

genes1 <- c("JAK1","IL15","CD36","STAT3","IRF3","TGFB1","TGFBR1","IFNGR1","TLR4","TLR2","MAPK14","PIAS1","PPARG","PPARA")

DotPlot(prot.combined, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltype.disease.region",
        idents = c("Monocytes_PN_PFC","DAM Microglia_PN_PFC","Homeostatic Microglia_PN_PFC",
                   "Monocytes_FTLD_PFC","DAM Microglia_FTLD_PFC","Homeostatic Microglia_FTLD_PFC"),
        scale = T,
        # split.by = "region_disease"
) + RotatedAxis() + coord_flip()


DotPlot(prot.combined, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltype.disease",
        idents = c("Homeostatic Microglia_PN","Homeostatic Microglia_FTLD",
                   "DAM Microglia_PN","DAM Microglia_FTLD",
                   "Monocytes_PN","Monocytes_FTLD"),
        scale = T,
        # split.by = "region_disease"
) + RotatedAxis() + coord_flip()

#####################           Volcanos
res <- read.xlsx("wilcox_MCX_FTLDxPN_Monocytes_DEGs.xlsx")
# res <- read.xlsx("wilcox_PFC_FTLDxPN_Homeostatic Microglia_DEGs.xlsx")
# res <- read.xlsx("wilcox_PFC_FTLDxPN_DAM Microglia_DEGs.xlsx")
# res <- read.xlsx("wilcox_PFC_FTLDxPN_CD8+ T cells_DEGs.xlsx")

colnames(res)[1] <- "genes"
# hav <- c(#"IFNG",
#   "FASLG","BATF","CD8A","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP","CX3CR1"
#   # "CD8B","IL7R","PRF1","IL2RB","GSDMB","KLRB1","LTB","TNF", "PDCD1","IRF2","CD3E",
# )
# 
# hav <- c("JAK1","IL15","STAT3","TGFB1","TGFBR1",#"CD36","IRF3",
#          "IFNGR2","TLR2","MAPK14",
#          "PIAS1","PPARG"
# )

# selected_geness <- c("IL15","IFNGR1","JAK1","JAK2","IRF1",  ### selection for DAM microglia by region and disease
#              "PIAS1","PPARA","PPARG","TLR4","TLR6",
#              "MAPK14","MAP2K1","MAP2K6","MAP3K3",  #"MAP2K4","MAP2K5","IRF2","MAP3K5","MAP4K3","MYD88","CD36",
#              "MAP4K5","MAPK8", ## "MAPK1","MAP4K4","MAP3K7","MAP3K14","MAP3K20","MAP3K4","MAP3K1",
#              "CREB1","CREBBP" #,"CD86","TLR2","IL15RA","IFNGR2","IRF3","STAT1","STAT2","STAT3","LRP1","SREBF1","SREBF2",
#              )

selected_geness <- c("IL15","JAK1","IRF3",  ### selection for monocytes by region and disease
            "PPARA","PPARG","CD36","TLR4","TLR6","SREBF2",
            "MAPK14","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K5",
            "MAP3K14","MAP3K20","MAP4K3","MAP4K5",  # "MAPK1","MAPK8","MAP4K4","MAP3K7","MAP3K4","MAP2K1",
            "CREB1","CREBBP"#,#,"CD86","IL15RA","LRP1","SREBF1","TLR2","PIAS1","IRF1", "IRF2","JAK2","STAT1","STAT2","STAT3","IFNGR1","IFNGR2","MYD88",
            )

# selected_geness <- c("IL15","JAK1","IRF2","IRF3",  ### selection for homeostatic microglia by region and disease
#            "PPARA","PPARG","TLR4","TLR6","SREBF2",
#            "MAPK14","MAP2K6","MAP3K3","MAP3K4", 
#            "MAP4K3","MAP4K5", # "MAPK1","MAPK8","MAP3K20","MAP4K4","MAP3K5","MAP3K7","MAP2K1","MYD88",
#            "CREB1","CREBBP"#,"CD86","IL15RA","LRP1","SREBF1","TLR2","IRF1","PIAS1","STAT1","STAT2","STAT3","IFNGR2","MAP3K14","MAP3K1","CD36","IFNGR1","JAK2",
#            )

# selected_geness <- c("STAT3","IRF2",#"IRF3", # signature for ----> monocytes
#                      "PIAS1","PPARG",#"MYD88",
#                      "MAP3K1","MAP3K3","MAP3K5","MAP4K5", # "MAPK1","MAPK8","MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP3K4","MAP2K6","MAP2K5","MAPK14","MAP2K1","MAP2K4",
#                      "CREB1",#,"CD86","IL15RA","IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","IRF1","PPARA","CD36","TLR2","TLR4","TLR6","CREBBP","LRP1","SREBF1","SREBF2",
#                      "IL15")

# selected_geness <- c("KLRB1",#"KLRK1","TNF","IFNG",   -------->> CD8 PFC
#                      "GSDMB","PDCD1",#"BATF","IL2RB", "IFNG-AS1", "IL2RG","GNLY",
#                      #"GZMA","GZMB","TYROBP","CD8B","CD8A", "NKG7","GZMH","PRF1","FASLG",
#                      "RUNX2","NCK2","CD2","RHOA","ITK","NFATC2","IKBKB","CD96","LAMP2", "FCMR","TOX","GZMK","LTB","IRF2","IL7R"
#                      )

# selected_geness <- c("KLRB1","TNF","IFNG",   ##-------->> CD8 MCX
#                      #"BATF","IL2RB", "IFNG-AS1", "IL2RG","GNLY",
#                      "GZMA","GZMB","TYROBP",#"CD8B","CD8A","NKG7","PDCD1","GZMH","PRF1","FASLG","KLRK1","GSDMB","GZMK","LTB","IKBKB","NFATC2","LAMP2", "TOX",
#                      "RUNX2","NCK2","CD2","RHOA","ITK","CD96","FCMR","IRF2","IL7R"
# )


res$neglog10p <- -log10(res$p_val)

# Assign color categories based on your thresholds
res$color <- "Not significant"
res$color[res$p_val < 0.05 & res$avg_log2FC > 0.05] <- "Upregulated"
res$color[res$p_val < 0.05 & res$avg_log2FC < -0.05] <- "Downregulated"

# Get coordinates for selected geness
label_data <- res[res$genes %in% selected_geness, ]

# Define vertical positions for stacked labels (right side)
x_label <- max(res$avg_log2FC) + 1.5
y_range <- range(res$neglog10p)
y_label <- seq(y_range[2], y_range[1], length.out = length(selected_geness))

# Data frame for label positions
stack_labels <- data.frame(
  genes = selected_geness,
  x_label = x_label,
  y_label = y_label
)

# Merge to get point coordinates for connector lines
stack_labels <- merge(stack_labels, label_data, by = "genes")

# Volcano plot with color by threshold and stacked labels
ggplot(res, aes(x = avg_log2FC, y = neglog10p)) +
  geom_point(aes(color = color), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Not significant" = "grey60", "Downregulated" = "lightblue")) +
  geom_segment(
    data = stack_labels,
    aes(x = avg_log2FC, y = neglog10p, xend = x_label, yend = y_label),
    color = "gray50", linewidth = 0.7, inherit.aes = FALSE
  ) +
  geom_label(
    data = stack_labels,
    aes(x = x_label, y = y_label, label = genes),
    hjust = 0, fontface = "bold", fill = "white", color = "black", label.size = 0.3, inherit.aes = FALSE
  ) +
  xlab(expression("Log"[2]*"(Fold-Change)")) +
  ylab(expression("-Log"[10]*"("*italic("p")*"-value)")) +
  ggtitle("Monocytes - FTLD x C (MCX)") +
  coord_cartesian(xlim = c(min(res$avg_log2FC) - 0.5, x_label + 2)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())




EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = selected_geness,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.05,
                # xlim = c(-7.5, 7.5),
                ylim = c(0, 30),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Monocytes - FTLD x C (MCX)",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                #maxoverlapsConnectors = NULL,
                #min.segment.length = 0.0000001,
                #directionConnectors = "x",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)















###################################        Myeloids analysis          ##################################
#################

myeloids <- subset(prot.combined, idents = c("5","12","46"))

### Normalize and scale merged obj
myeloids <- NormalizeData(myeloids, normalization.method = "LogNormalize")
myeloids <- FindVariableFeatures(myeloids, selection.method = "vst", nfeatures = 1000)
myeloids <- ScaleData(myeloids, vars.to.regress = c("percent.mt","nCount_RNA"), #latent.data = "nFeature_RNA", 
                      model.use = "linear")#, features = all.genes)
gc()

### Dimensionality reduction and integration
myeloids <- RunPCA(myeloids, npcs = 50)
gc()
ElbowPlot(myeloids, ndims = 50)
myeloids <- FindNeighbors(myeloids, dims = 1:30, reduction = "pca")
myeloids <- FindClusters(myeloids, resolution = 0.5, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
myeloids <- RunUMAP(myeloids, #umap.method = "umap-learn", 
                    dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(myeloids, reduction = "umap.unintegrated", raster = F,
        # ncol = 2,
        label = T,
        repel = T,
        #group.by = "patient"
        # split.by = "disease"
)#, combine = F)


## Save and Load data
saveRDS(object = myeloids, file = "myeloids_FTD_PFC_MCX_res0.5.Rds")
# #rm(myeloids)
myeloids <- readRDS("./myeloids_FTD_PFC_MCX_res0.5.Rds")

FeaturePlot(myeloids, raster = F, #pt.size = 1.5,
            features = "PRF1", 
            split.by = "region_disease",
            # ncol = 2,
            # min.cutoff = 1,
            # max.cutoff = 2.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")




## cell type annotation
myeloids[["RNA"]] <- JoinLayers(myeloids[["RNA"]])
Idents(myeloids) <- "seurat_clusters"
myeloids <- BuildClusterTree(myeloids, assay = "RNA", reduction = "umap.unintegrated", reorder = T)
prot.markers <- FindAllMarkers(myeloids, assay = "RNA", only.pos = T, min.pct = 0.35, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="myeloid.clusters.markers.xlsx")
prot.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5


DotPlot(myeloids, features = unique(top5$gene),
        col.max = 20,
        dot.scale = 10,
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


DotPlot(myeloids, features = c("dam2.markers1","dam1.markers1","Mono.markers1","MGND.markers1"),
        col.max = 20,
        dot.scale = 10,
        cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() #+ coord_flip()




DotPlot(myeloids, features = c("P2RY12","CSF1R","CD74","C3",#"CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",
                                    "IL1B","APOE","IRF8", #"CST7","CLEC7A","AIF1",
                                    "MS4A4A","CD163","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2","CD4"#,"GDA"#,"MARCO","SIGLEC1",
                                    #"CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7"#,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
),
col.max = 20,
dot.scale = 10, 
cluster.idents = T, #group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis() + coord_flip()



Idents(myeloids) <- "orig.ident"
mono.markers <- c("MS4A4A","CD163","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2"#,"CD4","GDA","MARCO","SIGLEC1"
)
DAM1 <- t(read.delim("DAM1.txt", header = F))
DAM2 <- t(read.delim("DAM2.txt", header = F))
MGND <- t(read.delim("MGND.txt", header = F))

myeloids <- AddModuleScore(
  myeloids,
  features = mono.markers,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Mono.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

myeloids <- AddModuleScore(
  myeloids,
  features = unique(DAM1),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "dam1.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

myeloids <- AddModuleScore(
  myeloids,
  features = unique(DAM2),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "dam2.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

myeloids <- AddModuleScore(
  myeloids,
  features = unique(MGND),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "MGND.markers",
  seed = 1,
  search = FALSE,
  slot = "data"
)

FeaturePlot(myeloids, features = "MGND.markers1", raster = F,# pt.size = 1.,
            # cols = "RdBu",
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 0.5,
            max.cutoff = 2.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")


monocytes <- subset(myeloids2, subset = Mono.markers1 > 0)

DotPlot(monocytes, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "clinical_dx",
        #idents = c("microglia"),
        scale = T,
        # split.by = "clinical_dx"
) + RotatedAxis() + coord_flip()


