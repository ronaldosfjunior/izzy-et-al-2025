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
#library(bluster)
library(pheatmap)
library(nichenetr)
library(tidyverse)
library(circlize)
library(ggpubr)


# setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/BPCells/")
# 
# file.dir <- "../h5_files/"
# files.set <- c(
#   "208_AD_4",
#   "221_AD_3",
#   "227_AD_2",
#   "228_AD_1",
#   "230_C_3",
#   "233_C_2",
#   "234_C_4",
#   "240_C_1"
# )
# data.list <- c()
# for (i in 1:length(files.set)) {
#   # path <- paste0(file.dir, sub("(.*_.*)_[0-9]", "\\1", files.set[i]), "_raw_feature_bc_matrix.h5")
#   # data <- open_matrix_10x_hdf5(path = path)
#   # write_matrix_dir(mat = data, dir = paste0(files.set[i], "_BP"))
#   ## Load in BP matrices
#   mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP"))
#   mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
#   dataset_name <- files.set[i]
#   condition <- sub(".*_(.*)_[0-9]", "\\1", files.set[i])
#   mat2 <- CreateSeuratObject(counts = mat, min.cells = 5, min.features = 500) %>% 
#     PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
#     subset(subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 20 & nCount_RNA < 20000) %>%
#   NormalizeData(normalization.method = "LogNormalize") %>%
#   FindVariableFeatures(selection.method = "vst")
#   #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)# %>%
#   #RunPCA(npcs = 50)
#   mat2$disease <- condition
#   mat2$patient <- dataset_name
#   mat2$cohort <- "Harvard"
#   mat2$cohort2 <- "MGH"
#   data.list[[i]] <- mat2
#   rm(mat)
#   rm(mat2)
#   #rm(data)
# }
# 
# setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/public_data/AD/PBMC/scRNAseq_PBMC_AD2/BPCells/")
# 
# ## Loop through h5 files and output BPCells matrices on-disk
# file.dir <- "../cellranger/"
# files.set <- c(
#   "GSM5145401_C1",
#   "GSM5145402_C2",
#   "GSM5145403_AD1",
#   "GSM5145404_AD2",
#   "GSM5145405_AD3",
#   "GSM5145406_AD4"
# )
# for (i in 1:length(files.set)) {
#   # path <- paste0(file.dir, files.set[i], "/outs/raw_feature_bc_matrix.h5")
#   # data <- open_matrix_10x_hdf5(path = path)
#   # write_matrix_dir(mat = data, dir = paste0(files.set[i], "_BP"))
#   ## Load in BP matrices
#   mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP"))
#   mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
#   dataset_name <- files.set[i]
#   condition <- sub(".*_(.*)[0-9]", "\\1", files.set[i])
#   mat2 <- CreateSeuratObject(counts = mat, min.cells = 5, min.features = 500) %>% 
#     PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
#     subset(subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 15 & nCount_RNA < 15000) %>%
#     NormalizeData(normalization.method = "LogNormalize") %>%
#     FindVariableFeatures(selection.method = "vst")
#   #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)
#   mat2$disease <- condition
#   mat2$patient <- dataset_name
#   mat2$cohort <- "Sichuan"
#   mat2$cohort2 <- "Sichuan"
#   data.list[[i+8]] <- mat2
#   rm(mat)
#   rm(mat2)
#   #rm(data)
# }
# 
# setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/yota/protollin/scRNAseq_clin_trial/Seurat/BPCells/PBMC/")
# file.dir <- "../../../cellranger/"
# 
# files.set <- c(
#   "GEM5",
#   "GEM7",
#   "GEM13",
#   "GEM15",
#   "GEM21",
#   "GEM22"
# )
# 
# for (i in 1:length(files.set)) {
#   ## Create binary matrices
#   # path <- paste0(file.dir, files.set[i], "/outs/filtered_feature_bc_matrix.h5")
#   # data <- open_matrix_10x_hdf5(path = path)
#   # write_matrix_dir(mat = data, dir = paste0(files.set[i], "_BP2"))
#   ## Load in BP matrices
#   mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP2"))
#   mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
#   dataset_name <- files.set[i]
#   mat2 <- CreateSeuratObject(counts = mat, min.cells = 5, min.features = 500) %>% PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#     subset(subset = nFeature_RNA < 5000 & percent.mt < 20 & nFeature_RNA > 500
#     ) %>% NormalizeData(normalization.method = "LogNormalize") %>%
#     FindVariableFeatures(selection.method = "vst")
#   mat2$patient <- dataset_name
#   mat2$disease <- "AD"
#   mat2$cohort <- "Harvard"
#   mat2$cohort2 <- "BWH"
#   data.list[[i+14]] <- mat2
#   rm(mat)
#   rm(mat2)
#   #rm(data)
# }
# 
# 
# setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/public_data/AD/PBMC/scRNAseq_PBMC_AD3/BPCells/")
# 
# ## Loop through h5 files and output BPCells matrices on-disk
# meta.data <- read.csv("meta.txt", header = T)
# file.dir <- "../cellranger/"
# files.set <- meta.data$GEM
# disease.state <- meta.data$Disease
# 
# for (i in 1:length(files.set)) {
#   # path <- paste0(file.dir, files.set[i], "/raw_feature_bc_matrix.h5")
#   # data <- open_matrix_10x_hdf5(path = path)
#   # write_matrix_dir(mat = data, dir = paste0(files.set[i], "_BP"))
#   ## Load in BP matrices
#   mat <- open_matrix_dir(dir = paste0(files.set[i], "_BP"))
#   mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
#   dataset_name <- files.set[i]
#   mat2 <- CreateSeuratObject(counts = mat, min.cells = 5, min.features = 500) %>% 
#     PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
#     subset(subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20 & nCount_RNA < 20000) %>%
#     NormalizeData(normalization.method = "LogNormalize") %>%
#     FindVariableFeatures(selection.method = "vst")
#   #SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)
#   mat2$disease <- disease.state[i]
#   mat2$patient <- dataset_name
#   mat2$cohort <- "Chicago"
#   mat2$cohort2 <- "Chicago"
#   data.list[[i+20]] <- mat2
#   rm(mat)
#   #rm(data)
# }
# 
# # Name layers
# files.set0 <- c(
#   "208_AD_4",
#   "221_AD_3",
#   "227_AD_2",
#   "228_AD_1",
#   "230_C_3",
#   "233_C_2",
#   "234_C_4",
#   "240_C_1",
#   "GSM5145401_C1",
#   "GSM5145402_C2",
#   "GSM5145403_AD1",
#   "GSM5145404_AD2",
#   "GSM5145405_AD3",
#   "GSM5145406_AD4",
#   "GEM5",
#   "GEM7",
#   "GEM13",
#   "GEM15",
#   "GEM21",
#   "GEM22"
# )
# files.set1 <- meta.data$GEM
# files.set <- c(files.set0, files.set1)
# names(data.list) <- files.set
# gc()

########  Merge everything
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/")

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

### Normalize and scale merged obj SCTransform
#prot.combined <- SCTransform(prot.combined, vars.to.regress = "percent.mt", return.only.var.genes = F)

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
prot.combined <- FindNeighbors(prot.combined, dims = 1:30, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 0.9, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
                              )
prot.combined <- RunUMAP(prot.combined, #umap.method = "umap-learn", 
                         dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

# DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
#         #ncol = 2,
#         label = T,
#         repel = T,
#         group.by = "seurat_clusters"
#         #split.by = "cohort2"
# )#, combine = F)


## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.Rds")
#rm(prot.combined)
#prot.combined <- readRDS("./obj_unintegrated_3000feats.Rds")


## Dimensionality reduction of integrated data
prot.combined <- RunHarmony(prot.combined, group.by.vars = "patient")
gc()
ElbowPlot(prot.combined, ndims = 50, reduction = "harmony")
prot.combined <- RunUMAP(prot.combined, dims = 1:40, reduction = "harmony", reduction.name = "umap")
prot.combined <- FindNeighbors(prot.combined, reduction = "harmony", dims = 1:40)
prot.combined <- FindClusters(prot.combined, resolution = 3, cluster.name = "harmony_clusters")
gc()

# # Using the bluster library, calculate the
# # 1-batch ARI
# batch_labels <- prot.combined@meta.data$cohort
# integrated_cluster_labels <- prot.combined@meta.data$seurat_clusters
# batch_ari <- 1 - pairwiseRand(batch_labels, integrated_cluster_labels,
#                               mode = "index", adjusted = TRUE)
# print(paste0("Batch mixing ARI: ", batch_ari))
# #> [1] "Batch mixing ARI: 0.952262112242732"

# ## Dimensionality reduction by UMAP
# prot.combined <- FindNeighbors(prot.combined, reduction = "umap", dims = 1:2)
# prot.combined <- FindClusters(prot.combined, resolution = 0.38, cluster.name = "umap_clusters")
# gc()

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F, 
        #ncol = 2,
        label = F,
        repel = T,
        group.by = "cohort",
        #split.by = "cohort"
) #, combine = F)

DimPlot(prot.combined, reduction = "umap", raster = F, 
        #ncol = 2,
        label = F,
        repel = T,
        #group.by = "seurat_clusters",
        #split.by = "age.range.disease2"
)#, combine = F)


## Save and Load data
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/")
saveRDS(object = prot.combined, file = "obj_harmony_patient_reg.out.mito.ncounts.v3c.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_harmony_patient_reg.out.mito.ncounts.v3c.Rds")
# saveRDS(object = prot.combined2, file = "no230_C_obj_harmony_patient_reg.out.mito.ncounts.v3b.Rds")
# #rm(prot.combined)
# prot.combined <- readRDS("./no230_C_obj_harmony_patient_reg.out.mito.ncounts.v3b.Rds")

meta.data <- read.csv("meta.txt", header = T)
prot.combined$age <- prot.combined$patient 
prot.combined$age.range <- prot.combined$patient
prot.combined$age.range2 <- prot.combined$patient
prot.combined$geno.disease <- prot.combined$patient
prot.combined$geno.disease2 <- prot.combined$patient
prot.combined$sex <- prot.combined$patient
prot.combined$geno <- prot.combined$patient

for (i in 1:length(meta.data$GEM)) {
  prot.combined$geno <- recode(prot.combined$geno, "meta.data$GEM[i] = meta.data$genotype[i]")
  prot.combined$sex <- recode(prot.combined$sex, "meta.data$GEM[i] = meta.data$sex[i]")
  prot.combined$age <- recode(prot.combined$age, "meta.data$GEM[i] = meta.data$Age[i]")
  prot.combined$age.range <- recode(prot.combined$age.range, "meta.data$GEM[i] = meta.data$age.range[i]")
  prot.combined$age.range2 <- recode(prot.combined$age.range2, "meta.data$GEM[i] = meta.data$age.range2[i]")
  prot.combined$geno.disease <- recode(prot.combined$geno.disease, "meta.data$GEM[i] = meta.data$geno.disease[i]")
  prot.combined$geno.disease2 <- recode(prot.combined$geno.disease2, "meta.data$GEM[i] = meta.data$geno.disease2[i]")
}
prot.combined$age.range.disease <- paste(prot.combined$age.range, prot.combined$disease, sep = "_")
prot.combined$age.range.disease2 <- paste(prot.combined$age.range2, prot.combined$disease, sep = "_")
prot.combined$sex.disease <- paste(prot.combined$sex, prot.combined$disease, sep = "_")

FeaturePlot(prot.combined, features = c("IFNG")
            , reduction = "umap",
            raster = F ,ncol = 2
            #,max.cutoff = 2.5
            , min.cutoff = 1.5
            #, blend = T
            , split.by = "age.range.disease2"
) +theme(legend.position = "right")

FeaturePlot(prot.combined, features = c("NKG7")
            , reduction = "umap",
            raster = F ,#ncol = 2
            #,max.cutoff = 2.5
            #, min.cutoff = 1.5
            #, blend = T
            #, split.by = "age.range.disease2"
) +theme(legend.position = "right")

irfs <- c("IRF1","IRF2",#"IRF3","IRF4","IRF5","IRF6",
          "IRF7","IRF8"#,"IRF9"
          )

RidgePlot(prot.combined, features = irfs, ncol = 2, sort = "decreasing", #group.by = "disease",
          idents = c("50_AD","47_AD","29_AD","7_AD","68_AD","71_AD","43_AD","69_AD","44_AD",
            "50_C","47_C","29_C","7_C","68_C","71_C","43_C","69_C","44_C"
                 
            ))
my_comparisons <- list(c("AD","C"))
my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
my_comparisons <- list(c("AD_E3/E3","C_E3/E3"),c("AD_E3/E4","C_E3/E4"),c("AD_E4/E4","C_E4/E4"),
                       c("AD_E4/E4","C_E3/E3"),c("AD_E4/E4","C_E3/E4"),c("AD_E3/E3","C_E4/E4")#,
                       #c("AD_E3/E3","C_E3/E3")
                       )
my_comparisons <- list(c("E2/E3","E3/E3"),c("E2/E3","E3/E4"),c("E2/E3","E4/E4"),
                       c("E3/E3","E3/E4"),c("E3/E3","E4/E4"),c("E3/E4","E4/E4"))

prot.combined$cohort2.disease <- paste(prot.combined$cohort2,prot.combined$disease,sep = "_")
my_comparisons <- list(c("Chicago_AD","Chicago_C"),c("MGH_AD","MGH_C"),c("Sichuan_AD","Sichuan_C"))
my_comparisons <- list(c("Chicago_AD","Chicago_C"),c("Harvard_AD","Harvard_C"),c("Sichuan_AD","Sichuan_C"))
MGH <- subset(prot.combined, subset = cohort2 == "MGH")
prot.combined2 <- subset(prot.combined, subset = patient != "230_C_3")
rm(prot.combined)
gc()

DotPlot(prot.combined,
        features = c("IL15","CD36"),
        group.by = "cohort.disease",
        idents = "Intermediate monocytes",
        col.max = 20,
        dot.scale = 10,
        scale = F
        )

p2 <- VlnPlot(prot.combined, features = "IL15",#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
        #split.by = "disease",
        group.by = "cohort.disease",
        pt.size = 0.05,
        raster = F,
        #ncol = 5,
        #slot = "counts",
        #add.noise = F,
        #log = T,
        #sort = "increasing",
        idents = ""#"Effector CD8+ T cells"
          #"Memory T CD8+"
        #"NKT-like CD8+"
          #"iNK"
          #"Classical monocytes"
          #c("1","2","20","41","21","67","62","65")
        #c("21","41","1","2","62","20","28","67","57") # NKs
          #c("80","18","48") # intermediate mono
           #c("19","58","81","78") # nonclassic.mono
        #c("50","47","29","7","68","71","43","69","44") #classic mono
) 
p2 <- p2 + scale_y_continuous(limits = c(0.000, max(p2[[1]][["data"]][["IL15"]])
                                         )) + stat_summary(fun = mean, geom = "point",size = 25, colour = "black", shape = 95) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p2$layers[[2]]$aes_params$alpha <- 0.1
p2

p2 <- VlnPlot(prot.combined, features = "IFNG",#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
              #split.by = "disease",
              group.by = "disease",
              pt.size = 0.001,
              raster = F,
              #ncol = 5,
              #slot = "counts",
              #add.noise = F,
              #log = T,
              #sort = "increasing",
              #idents = c("52","62") # IFNG producers
              #"Effector CD8+ T cells"
                #"Memory T CD8+"
                #"NKT-like CD8+"
                #"iNK"
                #"Classical monocytes"
                # c("1","2","20","41","21","67","62","65")
              #c("21","41","1","2","62","20","28","67","57") # NKs
              #c("80","18","48") # intermediate mono
              #c("19","58","81","78") # nonclassic.mono
              #c("50","47","29","7","68","71","43","69","44") #classic mono
) + theme(legend.position = "none") 
p2 <- p2 + scale_y_continuous(limits = c(0.000, max(p2[[1]][["data"]][["IFNG"]])+0.5
)) + stat_summary(fun = mean, geom = "point",size = 35, colour = "black", shape = 95)+ 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p2
# `21`= "iNK_1",
# `41`= "iNK_2",
# `67`= "iNK_3",
# `57`= "iNK_4",
# `62`= "iNK_5",
# `20`= "iNK_6",
# `1`= "mNK_1",
# `2`= "mNK_2",


################

genes1 <- c("JAK1","IL15","IL15RA","JAK2","CD36","CREB1","STAT1","STAT2","STAT3","IRF1","IRF2","IRF3","TGFB1","TGFBR2","TGFBR1",
            "TNF","IFNGR1","IFNGR2","TLR2","TLR4","TLR6","MAPK14","PIAS1","PPARG",
            "SPI1","CD14","JUN","JUNB","FOS","FOSL2","NFKB1","MYD88","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8") # Classical monocytes
genes1 <- c("PRF1","GZMA","GZMH","GZMK","GZMB","TOX","NKG7","IL7R","IRF2","IRF1","IRF3","IL2RG","IL2RB","STAT4","STAT1","STAT2","STAT3","JAK1","JAK2",
            "MAPK1","CD8A","CD8B","TRAC","TRBC1","CD28") # Effector CD8
genes1 <- c("CD86","CD80")

prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")
my_comparisons <- list(c("AD","C"))
my_comparisons <- list(c("50-70_male_AD","50-70_male_C"),c("71-90_male_AD","71-90_male_C"),
                       c("50-70_female_AD","50-70_female_C"),c("71-90_female_AD","71-90_female_C"))
my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "disease",
                #group.by = "sex.age2.disease",
                #group.by = "age.range.disease2",
                #group.by = "sex.disease",
                #idents = "Effector CD8+ T cells"
                idents = "Classical monocytes"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+0.5))+
    #scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+2))+
    #scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    #scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.disease.pdf"),
  #ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.age2.disease.pdf"),
  #ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.age.disease2.pdf"),
  #ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.disease.pdf"),
         plot = p1,
         width = 6, height = 8,
         #width = 12, height = 10,
         #width = 8, height = 8,
         #width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}


my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))
for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "sex.disease",
                #idents = "Classical monocytes"
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+ 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  #ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.disease.pdf"),
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.sex.disease.pdf"),
         plot = p1,
         width = 8,
         height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

######################

t.cell.cytotox <- c(
  "PRF1"
  ,"P2RX7"
  ,"PPP3CB"
  ,"B2M"
  ,"TAP2"
  ,"HLA-A"
  ,"RAB27A"
  ,"GZMM"
  ,"TRB"
  ,"TRA"
  ,"CRTAM"
  ,"HPRT1"
  ,"EMP2"
  ,"CTSC"
  ,"IL7R"
  ,"EBAG9"
  ,"CADM1"
  ,"CTSH"
  ,"MICA"
)

p1 <- VlnPlot(prot.combined, features = c(""),
        #split.by = "disease",
        group.by = "age.range.disease2",
        #group.by = "disease",
        pt.size = 0.05,
        raster = F,
        #ncol = 1,
        #slot = "counts",
        #add.noise = F,
        #log = T,
        #sort = "increasing",
        idents = "Effector CD8+ T cells"#c("50","47","29","7","68","71","43","69","44") #classic mono
        #c("80","18","48") # intermediate mono
        #c("19","58","81","78") # nonclassic.mono
        #c("62_Harvard_AD","62_Harvard_C","62_Chicago_AD","62_Chicago_C","62_Sichuan_AD","62_Sichuan_C")
) + stat_summary(fun = mean, geom = "point",size = 35, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, 8.5)) #+
  #geom_boxplot(colour = "grey")
p1$layers[[2]]$aes_params$alpha <- 0.2
p1

### celltypes records
library(ggstatsplot)
prot.combined$disease.patient <- paste(prot.combined$disease, prot.combined$patient, sep = "_")
num.cells <- as.data.frame(table(prot.combined$disease.patient, prot.combined$disease.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(prot.combined$disease.patient, prot.combined$seurat_clusters))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.clusters_disease_patient_prot.combined.csv")

df_freqs <- read.csv("freq.num.celltypes_disease_patient_prot.combined.csv",header = T)
plt <- ggbetweenstats(data = df_freqs,
                      x = geno.disease2,
                      y = Effector.CD8..T.cells,
                      p.adjust.method = "none",
                      type = "np"
)
plt
ggsave(filename = "Effector.CD8_freq.num_cells_np_disease_age2-vlnplot.pdf",
       plot = plt,
       width = 8,
       height = 8,
       device = "pdf")




df_freqs <- read.csv("freq.num.celltypes_disease_patient_prot.combined.csv", header = T)
celltypes2 <- colnames(df_freqs)[11:35]
# coluna1 <- c("AD_male","C_male","AD_female","C_female")
# coluna1 <- c("AD_50-70","C_50-70","AD_71-90","C_71-90")
coluna1 <- c("AD_50-70_male","C_50-70_male","AD_50-70_female","C_50-70_female",
             "AD_71-90_male","C_71-90_male","AD_71-90_female","C_71-90_female")

my_comparisons <- list(c("AD_50-70_male","C_50-70_male"),c("AD_50-70_female","C_50-70_female"),
                       c("AD_71-90_male","C_71-90_male"),c("AD_71-90_female","C_71-90_female"))
#my_comparisons <- list(c("AD_male","C_male"),c("AD_female","C_female"),c("AD_female","AD_male"),c("C_female","C_male"))
#my_comparisons <- list(c("AD_71-90","C_71-90"),c("AD_50-70","C_50-70"),c("AD_50-70","AD_71-90"),c("C_71-90","C_50-70"))
#my_comparisons <- list(c("AD","C"))


for (i in 1:length(celltypes2)){
  p1 <- ggplot(df_freqs, aes_string(x="factor(disease_sex_age2,levels=coluna1)", y=celltypes2[i], color="disease")) + 
    geom_violin(trim=T) + 
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=.3) + 
    geom_jitter(position=position_jitter(0.2)) +
    geom_boxplot(width=0.1) + 
    scale_color_brewer(palette="Dark2") + 
    #scale_fill_brewer(palette="Dark2") + 
    stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") +
    #stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
    stat_compare_means() +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") + # Add pairwise comparisons p-value
    theme_minimal() + rotate_x_text(angle = 45)
  #theme_classic()
  ggsave(filename = paste0("vln_stats_",celltypes2[i],"_by.disease_sex_age2.pdf"),
         plot = p1,
         width = 8,
         height = 6,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}
p1$layers[[2]]$aes_params$alpha <- 0.3



for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "sex.disease",
                #idents = "Classical monocytes"
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+ 
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  #ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.disease.pdf"),
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.sex.disease.pdf"),
         plot = p1,
         width = 8,
         height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

Idents(prot.combined) <- "celltypes"
markers <- FindAllMarkers(prot.combined, only.pos = TRUE, min.pct = 0.5) #%>%
  #group_by(cluster) #%>%
markers2 <- markers[markers$p_val_adj < 0.1,]
# dplyr::filter(avg_log2FC > 1)%>%
# slice_head(n = 50)
write.xlsx(as.data.frame(markers2), rowNames = T, file="wilcox_celltypes_all_markers_pct0.5_padj0.1.xlsx")

DotPlot(prot.combined, features = c("TNF","IFNG","KLRC1","NCAM1","IL2RB","IL7R","TBX21","EOMES",# iNKs
                                    "PRF1","GZMM","GZMH","GZMA","GZMB", # mNKs
                                    "LILRB1", "KLRB1", "ZBTB16", # NKT-like
                                    "CD3E","CD3D", # T cells
                                    "CD8A","CD8B","PTPRC","CCL5", #"CD244", # T cells
                                    "CD4","S100A4","SELL", # T cells
                                    "FOXP3", "IL2RA", # Tregs
                                    "TRDC","TRDV1","TRDV2","TRGC2","TRGV9","TRGC1", #gamma delta
                                    "ITGB2","PECAM1","IL3RA","LAMP1",#, #baso
                                    "CD64", "CD68", "CD71", "CCR5","ITGAM", # Macrophages
                                    "CD1C","ITGAX","FCER1A","CCR7","NRP1", # DC
                                    "CD14","FCGR3A", #monocytes
                                    "PF4", # platelets
                                    "CD19","IGKC","IGHM","CD27","CD1D","CD22","CD86","MS4A1","IGLC2","IGLC3","IGHD","CD79A","CD79B","AIM2", "BANK1","RALGPS2","TNFRSF13B", # B cells
                                    "IL4R","CXCR4", "BTG1", "TCL1A", "YBX3", # Naive B cells
                                    "COCH", "SSPN", "TEX9",  "TNFRSF13C", "LINC01781", # Memory B cells
                                    "LINC01857", # mature B cells
                                    "IGHA2","MZB1","TNFRSF17","DERL3","TXNDC5","POU2AF1","CPNE5","NT5DC2"# plasma cells
                                    # "IL10"#,"IL1B","IL15","IL7" #"IL2","IL3","IL6","IL4","IL22"
),
#cols = c("blue","blue"),#"blue"),#"red"),#"green","yellow","gray","pink","brown","lightblue"), 
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),#"Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"),
#c("Classical Mono_1","Classical Mono_2","Intermediate Mono","Nonclassical Mono", 
#"pDC_AD","pDC_C", 
#"mo-DC_AD","mo-DC_C"
#),
#idents = "CD8+ TEM",
# c("NK_4_C","NK_4_AD","NK_8_C","NK_8_AD","NK_21_C","NK_21_AD",
#  "CD8+ NKT-like_C", "CD8+ NKT-like_AD"#, "NK_AD", "NK_C"
# c("NK_AD","NK_C","Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C","Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"
#c("62","67","49","47"),
dot.scale = 10, 
cluster.idents = T, #group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis()

list1 <- c("NFKBIZ","CD83","BTG2","TNFRSF10B","RIPK2",#"NAMPT","PTGER4",
           "MAT2A","IL1B","ATF3","JUN","PLAUR")
list1 <- t(read.delim("NFKB_activated_targets.txt", header = F))
#list1 <- t(read.delim("NFKB_repressed_targets.txt", header = F))
DotPlot(prot.combined, features = c("MYD88","JUN","JUNB","CD14","IRF3"),#"FOS","FOSL2","NFKB1","IRAK1","RELA","TICAM1","TRAF6","TRAF3",),#list1,
col.max = 20, 
idents = "Classical monocytes",
dot.scale = 10,
cluster.idents = F, 
group.by = "condition",
#scale = F,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

DotPlot(prot.combined, features = c("CD8A","CD8B","IRF1","GZMA","PRF1","LTB","PDCD1","GZMH","GZMK","FASLG",#"IFNG","TNF",
                                    "TOX","IL2RB"#,"BATF","TYROBP","GSDMB"
),
#t.cell.cytotox,
#cols = c("blue","blue","blue"),
#cols = c("blue","red"),#"green","yellow","gray","pink","brown","lightblue"), 
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),#"Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"),
#c("Classical Mono_1","Classical Mono_2","Intermediate Mono","Nonclassical Mono", 
#"pDC_AD","pDC_C", 
#"mo-DC_AD","mo-DC_C"
#c("34","50","9","16","31","30","42"),
idents = "Effector CD8+ T cells",
# c("NK_4_C","NK_4_AD","NK_8_C","NK_8_AD","NK_21_C","NK_21_AD",
#  "CD8+ NKT-like_C", "CD8+ NKT-like_AD"#, "NK_AD", "NK_C"
# c("NK_AD","NK_C","Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C","Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"
# ),
dot.scale = 10,

cluster.idents = F, 
group.by = "age.range.disease2",
#scale = F,
#split.by = "disease"
) + RotatedAxis()

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "cluster.disease"

prot.combined$cell.age.disease <- paste(prot.combined$celltypes, prot.combined$age.range.disease2, sep = "_")
Idents(prot.combined) <- "cell.age.disease"


phagocytosis <- c(
  "FCN1",
  #"CD93",
  #"LYST",
  "NCF2",
  #"ANXA1",
  "ITGAM",
  #"RAB14",
  #"CLCN3",
  "PRKCD",
  "CORO1A",
  #"SPG11",
  "ITGAL",
  "P2RX7",
  "IRF8",
  #"PECAM1",
  #"TICAM2",
  "FCGR1A",
  #"CCR2",
  "NCF4",
  "CD14",
  "ICAM3",
  #"PAK1",
  #"ARHGAP25",
  "ITGB2",
  "VAV1",
  "BIN2",
  "RAC1",
  "HCK",
  #"CD36","NFKB1",
  #"ABL1"
  "SPI1","JUN","JUNB","FOS","FOSL2","MYD88"
)
phagocytosis2 <- c(
  "CD93",
  "LYST",
  "CLCN3",
  "CORO1A",
  "ITGAL",
  "P2RX7",
  "PECAM1",
  "FCGR1A",
  "CCR2",
  "ITGB2",
  "HCK",
  "CD36"
)
endocytosis <- c(
  #"CD93",
  #"LYST",
  "SNX17",
  #"CLCN3",
  "CORO1A",
  "ITGAL",
  "P2RX7",
  #"PECAM1",
  #"LIPA",
  #"DNAJC13",
  "FCGR1A",
  #"CAP1",
  #"CCR2",
  #"HEATR5A",
  "ASGR1",
  "PYCARD",
  #"DENND1A",
  #"FKBP15",
  "ITGB2",
  "HCK"#,
  #"LRRK2",
  #"CD36",
  #"RIN2",
  #"SNX10"
)
ROS <- c(
  "HVCN1",
  "ATP6V0D1",
  "NCF1",
  "ATP6V1G1",
  "ATP6V1B2",
  #"CYBB",
  "RAC2"
)

genes1 <- c("IL15","IL15RA","IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF1","IRF2","IRF3","PIAS1","PPARG","CD36",
            "TLR2","TLR4","TLR6","MAPK14","CREB1","TGFB1","TGFBR2","TGFBR1"#,
            #"TNF",
            #"SPI1","CD14","JUN","JUNB","FOS","FOSL2","NFKB1","MYD88","CD93"
            #"MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            #"MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8"
            ) # Classical monocytes

genes1 <- t(read.delim("NFKB_activated_targets_sel.txt", header = F))

my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)
DotPlot(prot.combined, features = genes1
#           c(#"CD8A","CD3E","IL2RB","GSDMB","LTB","PDCD1","FASLG","TNF","IL7R",
#           "IFNG","CD8B","BATF","KLRB1","IRF2","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP","IL2RG"#,"MAPK1","STAT4"
# )
,
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F,
# idents = c("Effector CD8+ T cells_71-90_AD","Effector CD8+ T cells_71-90_C",
#            "Effector CD8+ T cells_50-70_AD","Effector CD8+ T cells_50-70_C"),
idents = "Classical monocytes",
#scale = T,
group.by = "age.range.disease2"
#split.by = "disease"
) + RotatedAxis() + coord_flip()


stat1_targets_rep <- t(read.delim(file = "STAT1_targets_repression.txt",header = F))
stat1_targets <- t(read.delim(file = "STAT1_targets.txt",header = F))
stat1_regulators <- t(read.delim(file = "STAT1_activators.txt",header = F))
DotPlot(prot.combined, features = stat1_targets_rep,#c("TOX","CD8A","IRF1","GZMA","GZMH","GZMM","IFNG","TNF","PRF1","BATF","TYROBP","GSDMB"),#cd8act,
        cols = c("blue","blue"),
        #cols = c("blue","red"),#"green","yellow","gray","pink","brown","lightblue"), 
        col.max = 20, idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),#"Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"),
          #c("Classical Mono_1","Classical Mono_2","Intermediate Mono","Nonclassical Mono", 
          #"pDC_AD","pDC_C", 
          #"mo-DC_AD","mo-DC_C"
          #c("21","41","1","2","62","20","28","67","57") # NKs
          #c("80","18","48") # intermediate mono
          # c("19","58","81","78") # nonclassic.mono
          c("50","47","29","7","68","71","43","69","44"), #classic mono
          #c("36","63","52","9","16","31","67"),
        #idents = "CD8+ TEM",
        # c("NK_4_C","NK_4_AD","NK_8_C","NK_8_AD","NK_21_C","NK_21_AD",
        #  "CD8+ NKT-like_C", "CD8+ NKT-like_AD"#, "NK_AD", "NK_C"
        # c("NK_AD","NK_C","Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C","Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"
        # ),
        dot.scale = 10, 
        #cluster.idents = T, #group.by = "patho",
        #scale.by = "radius",
        scale = F,
        split.by = "disease"
) + RotatedAxis()


prot.combined <- RenameIdents(prot.combined, 
                              `50`= "Classical monocytes_1",
                              `47`= "Classical monocytes_2",
                              `29`= "Classical monocytes_3",
                              `7`= "Classical monocytes_4",
                              `68`= "Classical monocytes_5",
                              `71`= "Classical monocytes_6",
                              `43`= "Classical monocytes_7",
                              `69`= "Classical monocytes_8", # driven by two AD patients
                              `44`= "Classical monocytes_9",
                              `80`= "Intermediate monocytes_1", # driven by two AD patients, absent in HC
                              `18`= "Intermediate monocytes_2",
                              `48`= "Intermediate monocytes_3",
                              `36`= "Intermediate monocytes_4",
                              `19`= "Nonclassical monocytes_1",
                              `58`= "Nonclassical monocytes_2",
                              `78`= "Nonclassical monocytes_4", # driven by two AD patients, absent in HC
                              `21`= "iNK_1",
                              `41`= "iNK_2",
                              `67`= "iNK_3",
                              `57`= "iNK_4",
                              `62`= "iNK_5",
                              `20`= "iNK_6",
                              `1`= "mNK_1",
                              `2`= "mNK_2",
                              `22` = "Treg_1",
                              `61` = "Treg_2",
                              `33` = "Memory B cells_1",
                              `8` = "Naive B cells_1",
                              `60` = "Naive B cells_2",
                              `46` = "Naive B cells_3",
                              `63` = "Naive B cells_4",
                              `40` = "Naive B cells_5",
                              `38` = "Naive B cells_6",
                              `39` = "Cytotoxic B cells_1",
                              `37` = "mature B cells_1",
                              `70` = "mature B cells_2", # driven by two AD patients, absent in HC
                              `73` = "Plasma B cells_1",
                              `26` = "Cytotoxic B cells_1",
                              `23` = "Cytotoxic B cells_2", # looks like is driven by one or few patients
                              `65` = "Platelets_1",
                              `81` = "Platelets_2",
                              `77` = "Platelets_3",
                              `56` = "Platelets_4",
                              `45` = "mo-DC_1",
                              `53` = "pDC_1",
                              `76` = "pDC_2",
                              `75` = "pDC_3",
                              `3` = "NKT-like CD8+_1",
                              `30` = "NKT-like CD8+_2",
                              `52` = "NKT-like CD8+_3",
                              `72` = "NKT-like CD8+_4",
                              `12` = "NKT-like CD8+_5",
                              `28`= "NKT-like CD8+_6",
                              `17`= "NKT-like CD8+_7",
                              `25`= "NKT-like CD8+_8",
                              `64`= "NKT-like CD8+_9",
                              `34` = "NKT-like CD4+_1",
                              `66` = "NKT-like CD4+_2", # some AD patients has it increased
                              `82` = "NKT-like CD4+_3",
                              `31` = "Memory T CD8+_1",
                              `54` = "Memory T CD8+_2",
                              `5` = "Naive CD8+ T cells_1",
                              `74` = "Naive CD8+ T cells_2",
                              `9` = "Effector CD8+ T cells_1",
                              `14` = "Effector CD8+ T cells_2",
                              `79` = "Exhausted CD8+ T cells_1",
                              `59` = "Th1_1",
                              `32` = "Th1_2",
                              `51` = "Th1_3",
                              `49` = "Th1_4",
                              `16` = "gamma-delta T cells_1",
                              `24` = "Memory T CD4+_1",
                              `55` = "Memory T CD4+_2",
                              `0` = "Memory T CD4+_3",
                              `6` = "Memory T CD4+_4",
                              `10` = "Cytotoxic CD4+ T cells_1",
                              `15` = "Cytotoxic CD4+ T cells_2",
                              `4` = "Naive CD4+ T cells_1",
                              `11` = "Naive CD4+ T cells_2",
                              `42` = "Naive CD4+ T cells_3",
                              `27` = "Naive CD4+ T cells_4",
                              `35` = "Naive CD4+ T cells_5",
                              `13` = "Naive CD4+ T cells_6"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))


Idents(prot.combined) <- "seurat_clusters"
prot.combined <- RenameIdents(prot.combined, 
                              `50`= "Monocytes_1",
                              `47`= "Monocytes_2",
                              `29`= "Monocytes_3",
                              `7`= "Monocytes_4",
                              `68`= "Monocytes_5",
                              `71`= "Monocytes_6",
                              `43`= "Monocytes_7",
                              `69`= "Monocytes_8",
                              `44`= "Monocytes_9",
                              `80`= "Monocytes_1",
                              `18`= "Monocytes_2",
                              `48`= "Monocytes_3",
                              `36`= "Monocytes_4",
                              `19`= "Monocytes_1",
                              `58`= "Monocytes_2",
                              `78`= "Monocytes_4",
                              `21`= "NK_1",
                              `41`= "NK_2",
                              `67`= "NK_3",
                              `57`= "NK_4",
                              `62`= "NK_5",
                              `20`= "NK_6",
                              `1`= "NK_1",
                              `2`= "NK_2",
                              `22` = "CD4+ T cells_1",
                              `61` = "CD4+ T cells_2",
                              `33` = "B cells_1",
                              `8` = "B cells_1",
                              `60` = "B cells_2",
                              `46` = "B cells_3",
                              `63` = "B cells_4",
                              `40` = "B cells_5",
                              `38` = "B cells_6",
                              `39` = "B cells_1",
                              `37` = "B cells_1",
                              `70` = "B cells_2",
                              `73` = "B cells_1",
                              `26` = "B cells_1",
                              `23` = "B cells_2",
                              `65` = "Platelets_1",
                              `81` = "Platelets_2",
                              `77` = "Platelets_3",
                              `56` = "Platelets_4",
                              `45` = "DC_1",
                              `53` = "DC_1",
                              `76` = "DC_2", 
                              `75` = "DC_3", 
                              `3` = "NKT-like_1",
                              `30` = "NKT-like_2",
                              `52` = "NKT-like_3",
                              `72` = "NKT-like_4",
                              `12` = "NKT-like_5",
                              `28`= "NKT-like_6",
                              `17`= "NKT-like_7",
                              `25`= "NKT-like_8",
                              `64`= "NKT-like_9",
                              `34` = "NKT-like_1",
                              `66` = "NKT-like_2",
                              `82` = "NKT-like_3",
                              `31` = "CD8+ T cells_1",
                              `54` = "CD8+ T cells_2",
                              `5` = "CD8+ T cells_1",
                              `74` = "CD8+ T cells_2",
                              `9` = "CD8+ T cells_1",
                              `14` = "CD8+ T cells_2",
                              `79` = "CD8+ T cells_1",
                              `59` = "CD4+ T cells_1",
                              `32` = "CD4+ T cells_2",
                              `51` = "CD4+ T cells_3",
                              `49` = "CD4+ T cells_4",
                              `16` = "gamma-delta T cells_1",
                              `24` = "CD4+ T cells_1",
                              `55` = "CD4+ T cells_2",
                              `0` = "CD4+ T cells_3",
                              `6` = "CD4+ T cells_4",
                              `10` = "CD4+ T cells_1",
                              `15` = "CD4+ T cells_2",
                              `4` = "CD4+ T cells_1",
                              `11` = "CD4+ T cells_2",
                              `42` = "CD4+ T cells_3",
                              `27` = "CD4+ T cells_4",
                              `35` = "CD4+ T cells_5",
                              `13` = "CD4+ T cells_6"
)

prot.combined$celltypes3 <- sub("(.*)_.*", "\\1", Idents(prot.combined))



## find and visualize the top gene expression markers for each cluster
## Create downsampled object to make visualization either
Idents(prot.combined) <- "orig.ident"
object_subset <- subset(prot.combined, cells = Cells(prot.combined[["RNA"]]), downsample = 10000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "RNA"
Idents(object_subset) <- "celltypes"
object_subset <- BuildClusterTree(object_subset, assay = "RNA", reduction = "umap", reorder = T)

# markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE) %>%
#   group_by(cluster) #%>%
# # dplyr::filter(avg_log2FC > 1)%>%
# # slice_head(n = 50)
# write.xlsx(as.data.frame(markers), rowNames = T, file="wilcox_clus_all_markers_IPSI_50pcs_res3.xlsx")
markers <- FindAllMarkers(object_subset, assay = "RNA", only.pos = TRUE, min.pct = 0.8)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "RNA", features = top5$gene)
DoHeatmap(object_subset, assay = "RNA", features = top5$gene, size = 3.0) + 
  theme(axis.text = element_text(size = 5.5)) + NoLegend()

Idents(prot.combined) <- "celltypes"
DotPlot(object_subset, features = unique(top5$gene), cluster.idents = F,
        dot.scale = 10, 
        col.max = 20) + RotatedAxis() + coord_flip()








############  Diff exp analysis
# prot.combined$cohort.disease <- paste(prot.combined$cohort, prot.combined$disease, sep = "_")
prot.combined$cohort2.disease <- paste(prot.combined$cohort2, prot.combined$disease, sep = "_")
# prot.combined <- RenameIdents(prot.combined,`41` = "21",
#                               #`62` = "21",
#                               #`65` = "21",
#                               `21` = "21"
#                               )
Idents(prot.combined) <- "seurat_clusters"
prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
prot.combined$celltypes.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")
#zk.combined2$celltypes.condition <- paste(zk.combined2$customclassif, zk.combined2$condition, sep = "_")
prot.combined$cluster.condition <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
# prot.combined$cluster.condition <- paste(Idents(prot.combined), prot.combined$disease, sep = "_")
prot.combined$cluster.cohort <- paste(Idents(prot.combined), prot.combined$cohort, sep = "_")
prot.combined$cluster.cohort.condition <- paste(Idents(prot.combined), prot.combined$cohort, prot.combined$disease, sep = "_")
#zk.combined2$cluster.condition <- paste(zk.combined2$seurat_clusters, zk.combined2$condition, sep = "_")
#prot.combined$gender <- sub(".*_(.*)", "\\1", prot.combined$group)
#prot.combined$treatment <- sub("(.*)_.*", "\\1", prot.combined$group)
Idents(prot.combined) <- "seurat_clusters"
Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "celltypes.disease"
Idents(prot.combined) <- "cluster.condition"
Idents(prot.combined) <- "cluster.cohort.condition"
Idents(prot.combined) <- "cluster.cohort"
Idents(prot.combined) <- "customclassif"

prot.combined2$cluster.condition <- paste(Idents(prot.combined2), prot.combined2$disease, sep = "_")
Idents(prot.combined2) <- "cluster.condition"
nums <- as.character(0:82)

#######   MAST with regression comparing AD x C separating males and females

prot.combined$celltypes.disease.sex <- paste(prot.combined$celltypes, prot.combined$sex.disease, sep = "_")
Idents(prot.combined) <- "celltypes.disease.sex"

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_female_AD"),
                              ident.2 = paste0(celltypes[i],"_female_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("age.range2","cohort2"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_females_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_male_AD"),
                              ident.2 = paste0(celltypes[i],"_male_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("age.range2","cohort2"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_males_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

#######   MAST with regression comparing AD x C separating early and late

prot.combined$celltypes.disease.age2 <- paste(prot.combined$celltypes, prot.combined$age.range.disease2, sep = "_")
Idents(prot.combined) <- "celltypes.disease.age2"

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_50-70_AD"),
                              ident.2 = paste0(celltypes[i],"_50-70_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("sex","cohort2"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_50-70_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_71-90_AD"),
                              ident.2 = paste0(celltypes[i],"_71-90_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("sex","cohort2"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_71-90_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

##########   MAST with regression comparing AD x C

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("age.range2","sex","cohort"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 0.001,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              #recorrect_umi = TRUE
  )
  #write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("no230_MAST_lat.vars.disease.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

Idents(prot.combined) <- "cluster.condition"
nums <- as.character(81:82)
for (i in 1:length(nums)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(nums[i],"_AD"), 
                              ident.2 = paste0(nums[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.3,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("age.range2","sex","cohort"),
                              min.cells.feature = 3,
                              min.cells.group = 3,
                              pseudocount.use = 1,
                              mean.fxn = NULL,
                              fc.name = NULL,
                              base = 2,
                              densify = FALSE,
                              recorrect_umi = TRUE
  )
  zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_clusters_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_clus",nums[i],"_DEGs.xlsx"))
  rm(zk.response0)
}


Idents(prot.combined) <- "cluster.condition"

zk.response0 <- FindMarkers(prot.combined, ident.1 = "52_AD", 
                            ident.2 = "52_C",
                            slot = "data",
                            assay = "RNA",
                            features = NULL,
                            logfc.threshold = 0,
                            test.use = "MAST",
                            min.pct = 0.3,
                            min.diff.pct = -Inf,
                            verbose = TRUE,
                            only.pos = FALSE,
                            max.cells.per.ident = Inf,
                            random.seed = 1,
                            latent.vars = c("age.range2","sex","cohort"),
                            min.cells.feature = 3,
                            min.cells.group = 3,
                            pseudocount.use = 1,
                            mean.fxn = NULL,
                            fc.name = NULL,
                            base = 2,
                            densify = FALSE,
                            recorrect_umi = TRUE
)
zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
write.xlsx(as.data.frame(zk.response0), rowNames = T, file="MAST_clusters_lat.vars.age2.sex.cohort.Harvard_pct0.0_ADxHC_clus.52_DEGs.xlsx")
rm(zk.response0)

Idents(prot.combined) <- "celltypes.disease"

zk.response0 <- FindMarkers(prot.combined, ident.1 = "Intermediate monocytes_AD", 
                            ident.2 = "Intermediate monocytes_C",
                            slot = "data",
                            assay = "RNA",
                            features = NULL,
                            logfc.threshold = 0,
                            test.use = "MAST",
                            min.pct = 0.1,
                            min.diff.pct = -Inf,
                            verbose = TRUE,
                            only.pos = FALSE,
                            max.cells.per.ident = Inf,
                            random.seed = 1,
                            latent.vars = c("age.range2","sex","cohort"),
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
write.xlsx(as.data.frame(zk.response0), rowNames = T, file="MAST_clusters_lat.vars.age2.sex.cohort.Harvard_pct0.0_ADxHC_clus.52_DEGs.xlsx")
rm(zk.response0)





###### pseudo-bulk

# obj_samp <- subset(prot.combined, subset = celltypes == "Classical monocytes")
# 
# bulk <- AggregateExpression(obj_samp, assays = "RNA", slot = "counts", 
#                             return.seurat = T, group.by = c(#"celltypes",
#                                                             "disease","patient","age.range2","sex")
# )
# #bulk.nk.de <- FindMarkers(object = bulk, ident.1 = "CD8+ NKT-like cells", ident.2 = "CD4+ NKT-like cells", test.use = "DESeq2")
# bulk$celltype.disease <- paste(bulk$celltypes, bulk$disease, sep = "_")
# Idents(bulk) <- "celltype.disease"
# #bulk$patient <- sub(".*_(.*)_.*", "\\1", Cells(bulk))
# bulk$celltype.disease.patient <- paste(bulk$celltypes, bulk$disease,bulk$patient, sep = "_")
# Idents(bulk) <- "celltype.disease.patient"
# list1 <- grep("Classical monocytes_AD",bulk$celltype.disease.patient)
# list1 <- bulk$celltype.disease.patient[1:42]
# 
# list2 <- grep("Classical monocytes_C",bulk$celltype.disease.patient)
# list2 <- bulk$celltype.disease.patient[43:70]a
# 
# bulk2 <- bulk[grep("Classical monocytes",bulk[["RNA"]]@Dimnames[[2]])]
# 
# bulk.nk.de <- FindMarkers(object = bulk, ident.1 = list1,
#                           ident.2 = list2,
#                           test.use = "DESeq2")
# 
# Idents(bulk) <- "disease"
# bulk.nk.de <- FindMarkers(bulk, ident.1 = "AD",
#                           ident.2 = "C",
#                           slot = "data",
#                           assay = "RNA",
#                           features = NULL,
#                           logfc.threshold = 0,
#                           test.use = "negbinom",
#                           min.pct = 0.0,
#                           min.diff.pct = -Inf,
#                           verbose = TRUE,
#                           only.pos = FALSE,
#                           max.cells.per.ident = Inf,
#                           random.seed = 1,
#                           latent.vars = c("age.range2","sex"),
#                           min.cells.feature = 3,
#                           min.cells.group = 3,
#                           pseudocount.use = 0.001,
#                           mean.fxn = NULL,
#                           fc.name = NULL,
#                           base = 2,
#                           densify = FALSE,
#                           #recorrect_umi = TRUE
#                           )
# write.xlsx(as.data.frame(bulk.nk.de), rowNames = T, file="pseudobulk_negbinom_ADxHC_classic.mono_DEGs.xlsx")




########## DEGs records  -- lollipop plot

celltypes <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD8+",
  "NKT-like CD4+",
  "Platelets",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Effector CD8+ T cells",
  # "Exhausted CD8+ T cells",
  "Treg",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Th1",
  "Cytotoxic CD4+ T cells",
  "gamma-delta T cells",
  "Naive B cells",
  "Memory B cells",
  "Cytotoxic B cells",
  "mature B cells",
  "Plasma B cells"
)

data <- data.frame(matrix(ncol = 2,nrow = 0))

for (i in 1:length(celltypes)){
  list1 <- read.xlsx(paste0("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  data[i,1] <- celltypes[i]
  data[i,2] <- length(row.names(list1))
  data[i,3] <- length(list1$avg_log2FC[list1$avg_log2FC > 0])
  data[i,4] <- length(list1$avg_log2FC[list1$avg_log2FC < 0])
  rm(list1)
}

colnames(data) <- c("clusters","DEGs","upregulated","downregulated")

coluna1 <- data$clusters
coluna1 <- rev(coluna1) 
data$clusters <- factor(data$clusters, levels=coluna1)

# myCol <- randomColor(24, luminosity="light")#, hue = "green")
# 
# n <- 24
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# myCol <- paletteer_c("grDevices::Dynamic", 24) 
# myCol <- paletteer_c("grDevices::Earth", 24) 
myCol <- paletteer_c("grDevices::BluYl", 24) 

p1 <- ggplot(data, aes(x=factor(clusters,levels=coluna1), y=DEGs)) +
  geom_segment( aes(x=factor(clusters,levels=coluna1), xend=factor(clusters,levels=coluna1), y=0, yend=DEGs), color="skyblue") +
  geom_point( color=myCol,#data$color1,
              size=4, alpha=1) +
  theme_light() +
  scale_y_continuous(position = "right")+
  xlab("") +
  theme(axis.title.y = element_text(face="bold", color = "black", size=12),
        axis.text.x = element_text(face="bold", color = "black",size=11),
        axis.text.y = element_text(face="bold",color = "black",size=11)) +
  #coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

p1 <- p1 + coord_flip()
p1

####################################  Horizontal Bar plot ################################

axis_margin <- 5.5

p1 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x=upregulated)) +
  geom_col(fill = "darkorange") +
  scale_x_reverse() + 
  scale_y_discrete(position = "right") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
  ) + theme(
    ## Set background color to white
    panel.background = element_rect(fill = "white"),
    ## Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    ## Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    ## Remove the title for both axes
    #axis.title = element_blank(),
    ## Only left line of the vertical axis is painted in black
    #axis.line.y.left = element_line(color = "black"),
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
    ## But customize labels for the horizontal axis
    #axis.text.x = element_text(family = "Econ Sans Cnd", size = 10)
  )

p2 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x=downregulated)) +
  geom_col(fill = "skyblue") + scale_x_continuous(
    limits = c(0, 1500),
    breaks = seq(0, 1500, by = 500),
    #expand = c(0, 0), # The horizontal axis does not extend to either side
    #position = "bottom"  # Labels are located on the top
  ) +
  theme(
    axis.title.y = element_blank(),
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
  ) + theme(
    ## Set background color to white
    panel.background = element_rect(fill = "white"),
    ## Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    ## Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    ## Remove the title for both axes
    #axis.title = element_blank(),
    ## Only left line of the vertical axis is painted in black
    #axis.line.y.left = element_line(color = "black"),
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
    ## But customize labels for the horizontal axis
    #axis.text.x = element_text(family = "Econ Sans Cnd", size = 10)
  )

ggarrange(p1, p2)





####################################

NK41b <- subset(prot.combined, idents = c("41","62","21","65","2","1"))
NK41b <- subset(NK41b, subset = IFNG > 0.1)
NK41b$cluster.condition <- paste(Idents(NK41b), NK41b$disease, sep = "_")
Idents(NK41b) <- "cluster.condition"
NK41b[["RNA"]] <- JoinLayers(NK41b[["RNA"]])
zk.response0 <- FindMarkers(NK41b, ident.1 = "41_AD", 
                            ident.2 = "41_C",
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
write.xlsx(as.data.frame(zk.response0), rowNames = T, file="wilcox_ADxHC_NK.clus41_IFNG_DEGs.xlsx")
rm(zk.response0)






#####################      Save V5 as V3

# prot.combined[["RNA"]] <- JoinLayers(prot.combined[["RNA"]])
# prot.combined[["RNA"]]$scale.data <- NULL
# prot.combined[["RNA3"]] <- as(object = prot.combined[["RNA"]], Class = "Assay")
# prot.combined2 <- prot.combined
# 
# prot.combined2[["RNA"]] <- prot.combined2[["RNA3"]]
# prot.combined2[["RNA3"]] <- NULL
# 
# ## Save and Load data
# saveRDS(object = prot.combined, file = "obj_harmony_patient_reg.out.mito.ncounts.v3.Rds")
# #rm(prot.combined)
# prot.combined <- readRDS("./obj_harmony_patient_reg.out.mito.ncounts.v3.Rds")
# prot.combined2 <- ScaleData(prot.combined2)
# gc()
# Idents(prot.combined2) <- "cohort2"
# prot.combined2 <- subset(prot.combined2, idents = c("MGH","Sichuan","Chicago"))
# Idents(prot.combined2) <- "seurat_clusters"
# prot.combined2$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined2))
# prot.combined2$celltypes.disease <- paste(prot.combined2$celltypes, prot.combined2$disease, sep = "_")
# Idents(prot.combined2) <- "celltypes.disease"


# for (i in 1:length(celltypes)) {
#   zk.response0 <- FindMarkers(prot.combined2, ident.1 = paste0(celltypes[i],"_AD"), 
#                               ident.2 = paste0(celltypes[i],"_C"),
#                               slot = "data",
#                               assay = "RNA",
#                               features = NULL,
#                               logfc.threshold = 0,
#                               test.use = "wilcox",
#                               min.pct = 0.0,
#                               min.diff.pct = -Inf,
#                               verbose = TRUE,
#                               only.pos = FALSE,
#                               max.cells.per.ident = Inf,
#                               random.seed = 1,
#                               latent.vars = NULL,
#                               min.cells.feature = 3,
#                               min.cells.group = 3,
#                               pseudocount.use = 1,
#                               mean.fxn = NULL,
#                               fc.name = NULL,
#                               base = 2,
#                               densify = FALSE,
#                               recorrect_umi = TRUE
#   )
#   write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("wilcox_ADxHC_",celltypes[i],"_DEGs.xlsx"))
#   rm(zk.response0)
# }
# zk.response0 <- FindMarkers(prot.combined2, ident.1 = c("iNK_AD","mNK_AD"),#,"NKT-like CD8+_AD","NKT-like CD4+_AD"), 
#                             ident.2 = c("iNK_C","mNK_C"),#,"NKT-like CD8+_C","NKT-like CD4+_C"),
#                             slot = "data",
#                             assay = "RNA",
#                             features = NULL,
#                             logfc.threshold = 0,
#                             test.use = "wilcox",
#                             min.pct = 0.0,
#                             min.diff.pct = -Inf,
#                             verbose = TRUE,
#                             only.pos = FALSE,
#                             max.cells.per.ident = Inf,
#                             random.seed = 1,
#                             latent.vars = NULL,
#                             min.cells.feature = 3,
#                             min.cells.group = 3,
#                             pseudocount.use = 1,
#                             mean.fxn = NULL,
#                             fc.name = NULL,
#                             base = 2,
#                             densify = FALSE,
#                             recorrect_umi = TRUE
# )
# write.xlsx(as.data.frame(zk.response0), rowNames = T, file="wilcox_ADxHC_NKs_DEGs.xlsx")
# rm(zk.response0)

############## Identify cell types with scType

# lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); 
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# 
# # DB file
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
# tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# 
# # prepare gene sets
# gs_list = gene_sets_prepare(db_, tissue)
# 
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = prot.combined2[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
# 
# # merge by cluster
# cL_results = do.call("rbind", lapply(unique(prot.combined2@meta.data$seurat_clusters), function(cl){
#   es.max.cl = sort(rowSums(es.max[ ,rownames(prot.combined2@meta.data[prot.combined2@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
#   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(prot.combined2@meta.data$seurat_clusters==cl)), 10)
# }))
# sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])
# 
# prot.combined2@meta.data$customclassif = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   prot.combined2@meta.data$customclassif[prot.combined2@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# 
# # Visualization on UMAP
# DimPlot(prot.combined2, reduction = "umap", label = TRUE, repel = TRUE, 
#         group.by = 'customclassif', raster = F, 
#         #split.by = "condition"
# )   



#############################   Nichenet   #########################
library(nichenetr)
library(tidyverse)
library(circlize)

### subset by highly variable genes
#DefaultAssay(prot.combined) <- "RNA"
#prot.combined <- FindVariableFeatures(prot.combined)
#prot.combined2 <- subset(prot.combined, features = prot.combined@assays[["RNA"]]@var.features)

### start analysis
organism = "human"

#if(organism == "human"){
#  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
#} else if(organism == "mouse"){
#  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#}

## if connection does not work, add references by --->
setwd("/home/patrick/ref_data/nichenet_ref/human/")
lr_network = readRDS("lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/")

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

Idents(prot.combined) <- "seurat_clusters"

### circus plot attempt
## iNK ---> Classical Monocytes
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("21","41","62","20","67","57"#,"65"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks,
  top_n_ligands = 130)

## assign cluster
ligand_type_indication_df <- assign_ligands_to_celltype(prot.combined,
                                                        nichenet_output$top_ligands,
                                                        celltype_col = "celltypes") 

ligand_type_indication_df %>% head()
ligand_type_indication_df$ligand_type %>% table()

head(nichenet_output$ligand_target_df)
active_ligand_target_links_df <- nichenet_output$ligand_target_df
active_ligand_target_links_df$target_type <- "AD-DE" # needed for joining tables
circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
                                           active_ligand_target_links_df,
                                           cutoff = 0.40) 

head(circos_links)


ligand_colors <- c("General" = "#377EB8", "Exhausted CD8+ T cells" = "#4DAF4A", 
                   "Plasma B cells" = "#984EA3",
                   "Platelets" = "#FF7F00"#, 
                   #,"DC" = "#FFFF33"
                   , "pDC" = "#F781BF"
                   #,"CD8 T"= "#E41A1C"
                   ) 
target_colors <- c("AD-DE" = "#999999") 

vis_circos_obj <- prepare_circos_visualization(circos_links,
                                               ligand_colors = ligand_colors,
                                               target_colors = target_colors,
                                               celltype_order = NULL) 
make_circos_plot(vis_circos_obj, transparency = TRUE,  args.circos.text = list(cex = .5)) 

####### Analysis by cell types
Idents(prot.combined) <- "celltypes"

## NK ---> Monocytes
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("iNK","mNK"#,"NKT-like CD8+"#,"NKT-like CD4+"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## Monocytes ---> CD8
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("Naive CD8+ T cells"#,"Effector CD8+ T cells","Exhausted CD8+ T cells","Memory T CD8+"
  ), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## Monocytes ---> CD4
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("Memory T CD4+","Cytotoxic CD4+ T cells","Naive CD4+ T cells","Th1", "Treg"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks,
  top_n_ligands = 20)


## All ---> Classical Monocytes
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("Classical monocytes"),#"Intermediate monocytes","Nonclassical monocytes"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("iNK","mNK","NKT-like CD8+","NKT-like CD4+",
             "Memory T CD4+","Cytotoxic CD4+ T cells","Naive CD4+ T cells","Th1", "Treg",
             "Naive CD8+ T cells","Effector CD8+ T cells","Exhausted CD8+ T cells","Memory T CD8+",
             "Memory B cells","Naive B cells","Cytotoxic B cells","mature B cells","Plasma B cells","Cytotoxic B cells",
             "Platelets","mo-DC","pDC","gamma-delta T cells"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## Monocytes ---> NK NKT
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("iNK","mNK","NKT-like CD8+"#,"Memory T CD8+"
  ), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)





###### Analysis by clusters
## NK ---> Monocytes
## regular pipe clusters
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined2, 
  receiver = c("80","18","48","19","58","81","78","50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("67","57","28","21","41","62","1","2","20"#,"65"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


## NK ---> intermediate Monocytes
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined2, 
  receiver = c("80","18","48"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("21","41","1","2","62","20","28","67","57"#,"65"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## NK ---> Nonclassical Monocytes
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined2, 
  receiver = "19",#c("19","58","81","78"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("21","41","1","2","62","20","28","67","57"#,"65"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## NK ---> Classical Monocytes
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined2, 
  receiver = c("50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("21","41","1","2","62","20","28","67","57"#,"65"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


## NK ---> Monocytes
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined2, 
  receiver = "",#c("80","18","48","19","58","81","78","50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("67","57","28","21","41","65","1","2","20"#,"62"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


## NK IFNG ---> Monocytes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("50","47","29","7","68","71","43","69","44","19","58","78","80","18","48","36"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("41","52","62"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)

## NKT CD8 ---> Monocytes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("50","47","29","7","68","71","43","69","44","19","58","78","80","18","48","36"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("52"#,"3","30","72","12","28","17","25","64"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


## Monocytes ---> nk nkt
## regular pipe
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = all.cd8,#c("80","18","48","19","58","81","78","50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = classic,#c("67","57","28","21","41","65","1","2","20"),#,"62"
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)



## plots
nichenet_output$ligand_activities
write.xlsx(as.data.frame(nichenet_output$ligand_activities), rowNames = T, file="Monocytes_NK.NKT_nichenet_ligands.xlsx")
nichenet_output$top_ligands
nichenet_output$ligand_expression_dotplot + RotatedAxis()
nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_target_heatmap
nichenet_output$ligand_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  
                                                             high = "royalblue", breaks = c(0,0.0045,0.009)) + 
  xlab("AD response genes in NK") + ylab("Prioritized immmune cell ligands")
# DotPlot(prot.combined %>% subset(idents = c("10","23","38","30")), cols = c("royalblue","royalblue"),
#         features = nichenet_output$top_targets %>% rev(), split.by = "disease") + RotatedAxis()
nichenet_output$ligand_activity_target_heatmap
nichenet_output$ligand_receptor_heatmap


### cell types by clusters
eff.nkt.cd8 <- c("9","14","3","30","52","72","12","28","17","25","64")
all.cd8 <- c("9","14","31","54","5","74","79")
effect.cd8 <- c("9","14")
nktcd8 <- c("3","30","52","72","12","28","17","25","64")
NK <- c("1","2","67","57","21","41","62","20")
mNK <- c("1","2")
iNK <- c("67","57","21","41","62","20")
NK.NKT <- c("1","2","67","57","21","41","62","20","3","30","52","72","12","28","17","25","64")
classic <- c("50","47","29","7","68","71","43","44") # "69" has 1 cell in HC
nonclassic  <- c("19","58") # ,"78" absent in HC
intermediate  <- c("18","48","36") # "80", absent in HC
mono  <- c("50","47","29","7","68","71","43","44","19","58","18","48","36") # "69" has 1 cell in HC 80 and 78 absent in HC

####### Manual pipe
## receiver
receiver = all.cd8
expressed_genes_receiver = get_expressed_genes(receiver, prot.combined, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = classic
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, prot.combined, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
seurat_obj_receiver= subset(prot.combined, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["disease", drop=TRUE]])
condition_oi = "AD"
condition_reference = "C" 
DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% dplyr::filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% dplyr::mutate(rank = rank(dplyr::desc(aupr_corrected)))
#write.xlsx(as.data.frame(ligand_activities), rowNames = T, file="NK_mono.classical_nichenet_ligands.xlsx")
best_upstream_ligands = ligand_activities %>% top_n(121, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
best_upstream_ligands = file1$ligands
DotPlot(prot.combined, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()  ### Plot!!!!


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, 
                                                                 ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target2 <- data.frame(row.names = row.names(vis_ligand_target), targetnums = rowSums(vis_ligand_target != 0),
                                 mean.interaction = rowMeans(vis_ligand_target), median.interaction = rowMedians(vis_ligand_target)#,
                                 # mean.nozero = rowMeans(vis_ligand_target[vis_ligand_target != 0], na.rm = TRUE)
                                 )
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="all_classic_mono_2_all.CD8_nichenet_vis_ligand_target.xlsx")
# vis_ligand_target = vis_ligand_target[1:113,1:140]
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                                                                    color = "purple",legend_position = "top", y_axis = T,
                                                                    x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic", size = 4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.125,0.25))#c(0,0.0045,0.0090))
p_ligand_target_network ###  Plot!!!!

file1 <- read.xlsx("all_classic_mono_2_all.CD8_nichenet_vis_ligand_target2.xlsx")
colnames(file1)[1] <- "ligands"
colunas1 <- file1$ligands[15:1]
mid <- 0
ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(ligands, levels = colunas1), 
                                      size = targetnums, fill = mean.interaction), alpha = 1, shape = 21) +
  scale_size(range = c(3, 8), name = expression("Target genes"), breaks = c(20,60,100)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    #axis.text.x = element_text(face="bold", color = "black",size=12), 
    axis.text.y = element_text(face="bold",color = "black",size=10),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Classical onocytes -> CD8+ T cells") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$mean.interaction)), 
                       name = expression("Regulatory potential"#"Log"[2]*"(Fold-Enrichment)"
                       )) 


lr_network_top = lr_network %>% dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                                                                                        x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network ###  Plot!!!!


# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(prot.combined) %>% levels() %>% intersect(sender_celltypes) %>% 
  lapply(get_lfc_celltype, seurat_obj = prot.combined, condition_colname = "disease", 
         condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = NULL) %>% 
  reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
DE_table_all[is.na(DE_table_all)] = 0
# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% dplyr::rename(ligand = test_ligand) %>% left_join(DE_table_all %>% dplyr::rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0
# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]
colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", 
                                                                 low_color = "midnightblue",#mid_color = "white", mid = median(vis_ligand_lfc), 
                                                                 high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + 
  theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc  ###  Plot!!!!
# change colors a bit to make them more stand out
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue","lightblue","grey99","indianred1","firebrick1","red"),
                                                   #values = c(0,0.1,0.2,0.25, 0.40, 0.7,1),
                                                   #values = c(0,0.25,0.47,0.53,0.67,0.85,1) ## cd8 and NKs receiver from monocytes all main markers
                                                   values = c(0,0.38,0.55,0.67,0.78,0.9,1)
                                                   #limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1)
                                                   )
p_ligand_lfc


# ligand activity heatmap
ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()
vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% 
  make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                      x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_aupr ## Plot!!
# ligand expression Seurat dotplot
order_ligands_adapted <- str_replace_all(order_ligands, "\\.", "-")
rotated_dotplot = DotPlot(prot.combined %>% subset(seurat_clusters %in% sender_celltypes), # subset(seurat_clusters %in% sender_celltypes), 
                          features = order_ligands_adapted, 
                          cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), 
                                                                  legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
rotated_dotplot  ## Plot!!
figures_without_legend = cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), 
                          axis.text.y = element_text(face = "italic", size = 9), 
                          axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))
legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))
combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot


# `50`= "Classical monocytes_1",
# `47`= "Classical monocytes_2",
# `29`= "Classical monocytes_3",
# `7`= "Classical monocytes_4",
# `68`= "Classical monocytes_5",
# `71`= "Classical monocytes_6",
# `43`= "Classical monocytes_7",
# `69`= "Classical monocytes_8",
# `44`= "Classical monocytes_9",
# `80`= "Intermediate monocytes_1",
# `18`= "Intermediate monocytes_2",
# `48`= "Intermediate monocytes_3",
# `36`= "Intermediate monocytes_4",
# `19`= "Nonclassical monocytes_1",
# `58`= "Nonclassical monocytes_2",
# `78`= "Nonclassical monocytes_4",
# `21`= "iNK_1",
# `41`= "iNK_2",
# `67`= "iNK_3",
# `57`= "iNK_4",
# `62`= "iNK_5",
# `20`= "iNK_6",
# `1`= "mNK_1",
# `2`= "mNK_2",
# `3` = "NKT-like CD8+_1",
# `30` = "NKT-like CD8+_2",
# `52` = "NKT-like CD8+_3",
# `72` = "NKT-like CD8+_4",
# `12` = "NKT-like CD8+_5",
# `28`= "NKT-like CD8+_6",
# `17`= "NKT-like CD8+_7",
# `25`= "NKT-like CD8+_8",
# `64`= "NKT-like CD8+_9",










