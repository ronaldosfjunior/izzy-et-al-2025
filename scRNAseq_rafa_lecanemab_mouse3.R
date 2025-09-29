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
library(ComplexHeatmap)
library(igraph)
library(ggraph)
library(ggtext)


setwd("/media/patrick/JANDERSON/Bioinfo/weiner_lab/rafa/scRNAseq_lecanemab/Seurat/")


file.dir <- "../2025_04_Saef_scRNAseq_Somen/cellranger_outputs/"
data.list <- c()

files.set <- c("RR3_IC","RR1_IC_mLen", "C_Ctrl")

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i], "/filtered_feature_bc_matrix/")
  data <- Read10X(data.dir = path)
  dataset_name <- files.set[i]
  # condition <- sub("(.*)-.*", "\\1", files.set[i])
  # treatment1 <- treatment[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 500) %>% 
    PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20) %>%
    NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst")
  #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
  #RunPCA(npcs = 50)
  # mat$treatment <- treatment1
  mat$sample <- dataset_name
  data.list[[i]] <- mat
  rm(mat)
  rm(data)
}


setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/rafa/scRNAseq_LCMB/Seurat/")


file.dir <- "../cellranger/"

files.set <- c("M_Ctrl","RR7_IC","RR5_IC_mLen")

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i], "/filtered_feature_bc_matrix/")
  data <- Read10X(data.dir = path)
  dataset_name <- files.set[i]
  # condition <- sub("(.*)-.*", "\\1", files.set[i])
  # treatment1 <- treatment[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 500) %>% 
    PercentageFeatureSet(pattern = "^mt-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20) %>%
    NormalizeData(normalization.method = "LogNormalize") %>% FindVariableFeatures(selection.method = "vst")
  #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
  #RunPCA(npcs = 50)
  # mat$treatment <- treatment1
  mat$sample <- dataset_name
  data.list[[i+3]] <- mat
  rm(mat)
  rm(data)
}

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/rafa/scRNAseq_LCMB/Seurat_all/")

# Name layers
files.set <- c("RR3_IC","RR1_IC_mLen", "C_Ctrl","M_Ctrl","RR7_IC","RR5_IC_mLen")
names(data.list) <- files.set
# Merge layers and create seurat obj during merging 
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3500)
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
prot.combined <- FindVariableFeatures(prot.combined, selection.method = "vst", nfeatures = 3500)
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
prot.combined <- FindClusters(prot.combined, resolution = 0.9, cluster.name = "unintegrated_clusters")
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:40, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        # ncol = 2,
        label = T,
        repel = T,
        group.by = "seurat_clusters",
        split.by = "treatment"
)#, combine = F)

### add metadata
prot.combined$treatment <- prot.combined$sample 
sample1 <- c("RR3_IC","RR1_IC_mLen", "C_Ctrl","M_Ctrl","RR7_IC","RR5_IC_mLen")
sample2 <- c("IC","IC_mLen", "Ctrl","Ctrl","IC","IC_mLen")
for (i in 1:length(sample1)){
  prot.combined$treatment <- recode(prot.combined$treatment, "sample1[i] = sample2[i]")
}

prot.combined$treatment <- factor(prot.combined$treatment, levels = c("Ctrl","IC","IC_mLen"))

## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats.Rds")

prot.combined2 <- subset(prot.combined, subset = sample %in% c("RR3_IC","RR1_IC_mLen", "C_Ctrl"))
prot.combined2 <- subset(prot.combined2, subset = celltypes %in% c(
  "gdT",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "CD8+ NKT-like",
  "Th1",
  "iNK",
  "Memory CD8+ T cells",
  "Th2",
  "Exhausted CD8+ T cells",
  "Memory B cells",
  "Naive B cells",
  "NK"
))

## Dimensionality reduction of integrated data
prot.combined <- RunHarmony(prot.combined, group.by.vars = c("sample"))
gc()
ElbowPlot(prot.combined, ndims = 50, reduction = "harmony")
prot.combined <- RunUMAP(prot.combined, dims = 1:40, reduction = "harmony", reduction.name = "umap")
prot.combined <- FindNeighbors(prot.combined, reduction = "harmony", dims = 1:40)
prot.combined <- FindClusters(prot.combined, resolution = 0.9, cluster.name = "harmony_clusters")
gc()

DimPlot(prot.combined2, reduction = "umap", raster = F,
        # ncol = 2,
        label = T,
        repel = T,
        # group.by = "harmony_clusters",
        # split.by = "treatment"
) + theme(legend.position = "right")









## cell type annotation
prot.combined[["RNA"]] <- JoinLayers(prot.combined[["RNA"]])
prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap", reorder = T)
prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = T, min.pct = 0.35, logfc.threshold = 0.33) %>% group_by(cluster)
write.xlsx(as.data.frame(prot.markers), rowNames = T, file="all.celltypes.markers.harmony.xlsx")
prot.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top5


DotPlot(prot.combined, features = unique(top5$gene),
        col.max = 20,
        dot.scale = 10,
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() #+ coord_flip()

DotPlot(prot.combined, features = c(
  "Flt1", "Cldn5","Itm2a", #"Vtn","Vwf","Fam167b","Bmx","Clec1b",    #endothelial_cells
  # "Ambp","Higd1b","Cox4i2","Aoc3","Pde5a","Pth1r","P2ry14", "Abcc9","Kcnj8","Cd248", #Pericytes
  "P2ry12", "Csf1r","Cd74","C3","Hexb","C1qa","Aif1","Tmem119", "Il1b", #microglia "Cx3cr1","Cst3",
  "Ccr5","Itgam","Tfrc","Fcgr1", #macrophage 
  "Cx3cr1","Fcer1g","Lyz2","Ccr2","Ly6c1","Cd68","Cd14","Fcgr3","Il15","F13a1", #monocytes "Itgam",
  "Itgax","Cd1d1","Nrp1","Il3ra","Pecam1", # DC  "Itgb2","Lamp1","Bst2",
  "S100a9","Cxcr2","Mmp9", "S100a8","Il1r2","Retnlg", # Neutrophils 
  "Pf4", # platelets
  "Ifng","Tnf", # Th1
  "Gata3", # Th2
  "Il17a",
  "Tox", "Il2rg",  # iNKs 
  "Klrb1c","Il2rb", "Ncam1","Klrc1",
  "Prf1","Gzma","Gzmm","Gzmk","Gzme","Gzmb", # mNKs 
  "Il7r","Cd8a","Cd8b1",#"Ptprc",#"CCL5",#"CD244",
  "Cd4","S100a4","Sell", # T cells
  "Cd3e","Cd3d", "Ccr7",# T cells
  "Foxp3", "Il2ra", # Tregs
  "Trdc","Trgv2","Trgv4","Trgv5","Trdv4", #"Trgv7","Trdv1", #gamma delta
  "Cd19","Cd22","Ighg1","Igkc","Ighv1-26","Cd27","Cd86","Ms4a1" #,"IGLC2","IGLC3"  "Ighm",#, "IGLL5" b cells
),
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),
dot.scale = 10, 
cluster.idents = T, #group.by = "seurat_clusters",
#scale = F,
#split.by = "condition"
) + RotatedAxis() #+ coord_flip()

DotPlot(prot.combined, features = c(
  "Apoe","Gfap", "Aqp4", "Lcn2", "Gja1", "Slc1a2", "Fgfr3", "Nkain4",   #astrocytes
  "Slc17a6",  "Slc17a7",  "Nrgn", "Camk2a", "Satb2", "Col5a1", "Sdk2", "Nefm",    #excitatory_neurons
  "Slc32a1",  "Gad1", "Gad2", "Tac1", "Penk", "Sst",  "Npy",  "Mybpc1", "Pvalb", "Gabbr2",   #inhibitory_neurons
  "Flt1", "Cldn5", "Vtn", "Itm2a", "Vwf", "Fam167b", "Bmx", "Clec1b",    #endothelial_cells
  "Ambp",  "Higd1b", "Cox4i2", "Aoc3", "Pde5a",  "Pth1r",  "P2ry14", "Abcc9", "Kcnj8", "Cd248", #Pericytes
  "Olig2",  "Mbp",  "Mobp", "Plp1", "Mog",  "Cldn11", "Myrf", "Galc", "Ermn", "Mag",   #oligodendrocytes
  "Vcan", "Cspg4", "Pdgfra", "Sox10", "Neu4", "Pcdg15", "Gpr37l1", "C1ql1", "Cdo1", "Epn2",   #oligodendrocyte_precursor_cells
  "Iba1", "P2ry12", "Csf1r",  "Cd74", "C3", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia
  "Ccr2","Cd68","Cd11b","Cd14","Fcgr3" #monocytes
),
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),
dot.scale = 10, 
cluster.idents = T, #group.by = "seurat_clusters",
#scale = F,
#split.by = "condition"
) + RotatedAxis() #+ coord_flip()

DotPlot(prot.combined, features = c("P2ry12","Csf1r","Cd74","C3","Cst3","Hexb", "C1qa", "Cx3cr1","Tmem119", #"SLC2A5",
                                    "Il1b","Apoe","Irf8", "Cst7","Clec7a","Aif1",
                                    "Marco","Ms4a4a","Cd63","Siglec1","F13a1","Mafb","Tnfaip2","Il15","Asah1","Gas7","Pla2g7","Plxnd1","Emilin2","Cd4","Gda" #,
                                    # "CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7"#,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
),
col.max = 20,
dot.scale = 10, 
cluster.idents = T, #group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis() #+ coord_flip()


prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Endothelial cells",
                              `1`= "Naive CD8+ T cells",
                              `2`= "Monocytes",
                              `3`= "Homeostatic Microglia", 
                              `4`= "Th1",
                              `5`= "Th2",
                              `6`= "Memory CD8+ T cells",
                              `7`= "Neuroendocrine cells",
                              `8`= "Monocytes",
                              `9`= "CD8+ NKT-like",
                              `10`= "gdT",
                              `11`= "Memory CD8+ T cells",
                              `12`= "BM Macrophage",
                              `13`= "Exhausted CD8+ T cells",
                              `14`= "Monocytes",
                              `15`= "NK_1", 
                              `16`= "Effector CD8+ T cells",
                              `17`= "Homeostatic Microglia",
                              `18`= "Progenitor cells", 
                              `19`= "Oligodendrocytes",
                              `20`= "Naive B cells",
                              `21`= "Pericytes",
                              `22`= "Macrophages",
                              `23`= "MgND Microglia",
                              `24`= "Monocytes",
                              `25`= "Macrophages", 
                              `26`= "Naive CD8+ T cells",
                              `27`= "Endothelial cells",
                              `28`= "Th2",
                              `29`= "cDC",
                              `30`= "Memory B cells",
                              `31`= "iNK",
                              `32`= "Astrocytes",
                              `33`= "Exhausted CD8+ T cells",         
                              `34`= "Mesenchymal stem cells",
                              `35`= "Oligodendrocytes",
                              `36`= "Neutrophils",
                              `37`= "Endothelial cells",
                              `38`= "Endothelial cells",
                              `39` = "VSMC",
                              `40` = "gdT",
                              `41` = "iNK",
                              `42` = "Homeostatic Microglia"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"

## Save and Load data
saveRDS(object = prot.combined, file = "obj_integrated_3500feats_annotated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_integrated_3500feats_annotated.Rds")


cell.types <- c(
  "gdT",
  "Th2",
  "Th1",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "Memory CD8+ T cells",
  "CD8+ NKT-like",
  "NK",
  "iNK",
  "Exhausted CD8+ T cells",
  "Memory B cells",
  "Naive B cells"
)

cell.types2 <- c(
  "gdT",
  "gdT",
  "Th1",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "Memory CD8+ T cells",
  "CD8+ NKT-like",
  "NK",
  "progenitor lymphoid cell",
  "Exhausted CD8+ T cells",
  "Memory B cells",
  "Naive B cells"
)

prot.combined2$celltypes2 <- prot.combined2$celltypes 
for (i in 1:length(cell.types2)){
  prot.combined2$celltypes2 <- recode(prot.combined2$celltypes2, "cell.types[i] = cell.types2[i]")
}

cell.types <- c(
  "gdT",
  "Th1",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "Memory CD8+ T cells",
  "CD8+ NKT-like",
  "NK",
  "progenitor lymphoid cell",
  "Exhausted CD8+ T cells",
  "Memory B cells",
  "Naive B cells"
)

prot.combined2$celltypes2 <- factor(prot.combined2$celltypes2, levels = cell.types)


DotPlot(prot.combined2, features = c(
  #"Flt1", "Cldn5","Itm2a", #"Vtn","Vwf","Fam167b","Bmx","Clec1b",    #endothelial_cells
  # "Ambp","Higd1b","Cox4i2","Aoc3","Pde5a","Pth1r","P2ry14", "Abcc9","Kcnj8","Cd248", #Pericytes
  #"P2ry12", "Csf1r","Cd74","C3","Hexb","C1qa","Aif1","Tmem119", "Il1b", #microglia "Cx3cr1","Cst3",
  #"Ccr5","Itgam","Tfrc","Fcgr1", #macrophage 
  #"Cx3cr1","Fcer1g","Lyz2","Ccr2","Ly6c1","Cd68","Cd14","Fcgr3","Il15","F13a1", #monocytes "Itgam",
  #"Itgax","Cd1d1","Nrp1","Il3ra","Pecam1", # DC  "Itgb2","Lamp1","Bst2",
  #"S100a9","Cxcr2","Mmp9", "S100a8","Il1r2","Retnlg", # Neutrophils 
  #"Pf4", # platelets
  "Trdc","Trgv2","Trgv4","Trgv5","Trdv4","Cd4",
  "Ifng", "Cd3e","Cd3d","Cd8a","Cd8b1", #"Tnf", # Th1
  #"Gata3", # Th2
  #"Il17a",
   #"Il2rg",  # iNKs "Ncam1","Il2rb",
  "Klrb1c", "Klrc1",
  "Prf1","Gzma","Tox", #"Gzmm","Gzmk","Gzme","Gzmb", # mNKs 
  #"Ptprc",#"CCL5",#"CD244","Il7r",
   #"S100a4","Sell", # T cells
  # 
  "Ccr7",# T cells
  # "Foxp3", "Il2ra", # Tregs
   #"Trgv7","Trdv1", #gamma delta
  # "Cd19","Cd22","Cd27",
  "Cd86","Ms4a1","Igkc","Ighv1-26" #,"IGLC2","IGLC3"  "Ighm","Ighg1",#, "IGLL5" b cells
),
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),
dot.scale = 10, 
cluster.idents = F, group.by = "celltypes2",
#scale = F,
#split.by = "condition"
) + RotatedAxis() #+ coord_flip()


prot.combined$celltypes <- factor(prot.combined$celltypes, levels = cell.types)


FeaturePlot(prot.combined2, features = "Gzmb", raster = T, pt.size = 2.,
            cols = c("gray90","#9E1021"),
            split.by = "treatment",
            # ncol = 2,
            # min.cutoff = 2.,
            # max.cutoff = 3.0,
            reduction = "umap") + theme(legend.position = "right")

DimPlot(prot.combined2, reduction = "umap", raster = F,
        # ncol = 2,
        label = T,
        repel = T,
        group.by = "celltypes2"
        # split.by = "treatment"
)#, combine = F)

prot.combined$treatment <- factor(prot.combined$treatment, levels = c("Ctrl","IC","IC_mLen"))
my_comparisons <- list(c("IC","IC_mLen"),c("Ctrl","IC"),c("Ctrl","IC_mLen"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,3)]

gene <- "Il15"
p1 <- VlnPlot(prot.combined, features = gene,
              pt.size = 0.05, raster = F, group.by = "treatment", cols = colors,
              idents = c("")
              # idents = c("Effector CD8+ T cells")
              # idents = c("mNK")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000005, max(p1[[1]][["data"]][[gene]])+2.))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1



DotPlot(prot.combined, features = c(#"Ifng","Tnf", # "Ncam1","Il7r",
                                    #"Tox", "Il2rg",
                                    "Il2rb",   #  
                                    # "Klrb1c","Klrc1","Klrb1b",
                                    "Prf1","Gzma","Gzmm","Gzmk","Gzme","Gzmb", # mNKs 
                                    "Cd8a","Cd8b1",
                                    "Stat3","Ifngr1","Il18rap","Irf8" # ,"Klri2","Tnfaip3",
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "treatment",
idents = c("Effector CD8+ T cells","CD8+ NKT-like",
           "NK"),
scale = T,
# split.by = "treatment"
) + RotatedAxis() + coord_flip()


DotPlot(prot.combined, features = c("Il15","Il15ra","Irf2","Pias1","Tnf","Mapk14","Mapk1","Spi1","Pparg","Creb1","Stat1","Cd36","Tlr2","Tlr4","Jak1",
                                    "Jak2","Ifngr1","Ifngr2"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "treatment",
idents = c(
  "Monocytes" #,
  # "MgND Microglia" #,
  # "BM Macrophage" #,
  # "Macrophages" #,
  # "Homeostatic Microglia"
  ),
scale = T,
# split.by = "treatment"
) + RotatedAxis() + coord_flip()



IL15.sig <- c("Il15","Il15ra","Irf2","Pias1","Tnf","Mapk14","Mapk1","Spi1","Pparg" ,   "Creb1","Stat1","Cd36","Tlr2","Tlr4","Jak1",
              "Jak2","Ifngr1","Ifngr2"
) # classic Monocytes

prot.combined <- AddModuleScore(
  prot.combined,
  features = IL15.sig,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "IL15.signature",
  seed = 1,
  search = FALSE,
  slot = "data"
)

gene1 <- "IL15.signature1"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "treatment", cols = colors,
              idents = c("Endothelial cells"
                        #"Monocytes",
                        #"MgND Microglia",
                        #  "BM Macrophage",
                        # "Macrophages",
                        # "Homeostatic Microglia"
              )
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[gene1]])+0.5*max(p1[[1]][["data"]][[gene1]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.2
p1


cytotox <- c(#"Ifng",#"Tnf", # "Ncam1","Il7r",
             # "Tox", "Il2rg","Il2rb",   #  
             # "Klrb1c","Klrc1",
             # "Prf1","Gzma","Gzmm","Gzme","Gzmb", # mNKs  "Gzmk",
             # "Cd8a","Cd8b1"
  #"Ifng","Tnf", # "Ncam1","Il7r",
  #"Tox", "Il2rg",
  "Il2rb",   #  
  # "Klrb1c","Klrc1","Klrb1b",
  "Prf1","Gzma","Gzmm","Gzmk","Gzme","Gzmb", # mNKs 
  "Cd8a","Cd8b1",
  "Stat3","Ifngr1","Il18rap","Irf8" # ,"Klri2","Tnfaip3",
) # CD8

prot.combined <- AddModuleScore(
  prot.combined,
  features = cytotox,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Cytotoxicity",
  seed = 1,
  search = FALSE,
  slot = "data"
)

gene1 <- "Cytotoxicity1"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "treatment", cols = colors,
              idents = c(#"Effector CD8+ T cells",
                         "CD8+ NKT-like" #, 
                         #"NK"
                         )
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.00000, max(p1[[1]][["data"]][[gene1]])+0.35*max(p1[[1]][["data"]][[gene1]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.2
p1


############### fig UMAP color


cell.types <- c(
  "progenitor lymphoid cell",
  "gdT",
  "Th1",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "Memory CD8+ T cells",
  "CD8+ NKT-like",
  "NK",
  "Exhausted CD8+ T cells",
  "Memory B cells",
  "Naive B cells"
)

prot.combined2$celltypes2 <- factor(prot.combined2$celltypes2, levels = cell.types)

col1 <- c(#paletteer_c("ggthemes::Orange Light", 1), 
          paletteer_c("ggthemes::Gold-Purple Diverging", 1),
          paletteer_c("ggthemes::Blue", 2),
          paletteer_c("ggthemes::Orange-Gold", 6),
          paletteer_c("ggthemes::Purple", 2))

DimPlot(prot.combined2, reduction = "umap", raster = T, pt.size = 1.5, 
        cols = col1,
        label = T,
        repel = T,
        group.by = "celltypes2",
        # split.by = "cohort"
) + theme(legend.position = "right")

















