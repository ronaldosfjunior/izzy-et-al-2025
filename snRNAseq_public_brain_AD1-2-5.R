# calling packages
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(MAST)
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
library(openxlsx)
library(BiocParallel)
register(MulticoreParam(12))
library(harmony)
library(reticulate)
library(ggthemes)
options(future.globals.maxSize = 8e+09)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr) 
#library(bluster)
#library(pheatmap)

# ##################### AD brain 1 - Ip
# setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/scRNAseq_Ip/Seurat_new/")
# 
# file.dir <- "../cellranger_count/"
# files.set <- c(
#   "AD/GSM4775561_AD1/",
#   "AD/GSM4775562_AD2/",
#   "AD/GSM4775563_AD4/",
#   "AD/GSM4775564_AD5/",
#   "AD/GSM4775565_AD6/",
#   "AD/GSM4775566_AD8/",
#   "AD/GSM4775567_AD9/",
#   "AD/GSM4775568_AD10/",
#   "AD/GSM4775569_AD13/",
#   "AD/GSM4775570_AD19/",
#   "AD/GSM4775571_AD20/",
#   "AD/GSM4775572_AD21/",
#   "C/GSM4775573_NC3/",
#   "C/GSM4775574_NC7/",
#   "C/GSM4775575_NC11/",
#   "C/GSM4775576_NC12/",
#   "C/GSM4775577_NC14/",
#   "C/GSM4775578_NC15/",
#   "C/GSM4775579_NC16/",
#   "C/GSM4775580_NC17/",
#   "C/GSM4775581_NC18/"
# )
# 
# data.list <- c()
# 
# for (i in 1:length(files.set)) {
#   path <- paste0(file.dir, files.set[i])
#   data <- Read10X(data.dir = path)
#   dataset_name <- sub(".*[A-Z]/(GSM.*)/", "\\1", files.set[i])
#   mat <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200) %>% 
#     PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
#     subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5) %>%
#     NormalizeData(normalization.method = "LogNormalize") %>%
#     FindVariableFeatures(selection.method = "vst") %>% ScaleData()
#   mat$disease <- sub("(.*[A-Z])/GSM.*/", "\\1", files.set[i])
#   mat$patient <- dataset_name
#   mat$cohort <- "Hong Kong University"
#   data.list[[i]] <- mat
#   rm(mat)
#   rm(data)
# }
# 
# # # Name layers
# # names(data.list) <- files.set
# 
# # # Merge layers and create seurat obj during merging 
# # features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
# # prot.combined1 <- merge(data.list[[1]], y = data.list[2:length(data.list)], 
# #                        add.cell.ids = files.set, merge.data = T)
# # 
# # VariableFeatures(prot.combined1) <- features
# # rm(data.list)
# gc()
# 
# #####################  AD brain 2 - Tsai
# 
# setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/scRNAseq_tsai/Seurat_new/")
# data <- Read10X(data.dir = "../cellranger_processed/all/")
# 
# prot.combined2 <- CreateSeuratObject(counts = data,  min.cells = 3, min.features = 200) %>%
#   PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#   subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5) %>%
#   NormalizeData(normalization.method = "LogNormalize") %>%
#   FindVariableFeatures(selection.method = "vst") %>% ScaleData()
# 
# rm(data)
# gc()
# 
# ## add meta
# prot.combined2$patient <- sub(".*-(.*)", "\\1", prot.combined2@assays[["RNA"]]@data@Dimnames[[2]])
# prot.combined2$patients <- as.numeric(prot.combined2$patient)
# no_patho <- "C"
# early_patho <- "Early AD"
# late_patho <- "Late AD"
# patho <- "AD"
# no_patho.list <- scan("list1.txt")
# early_patho.list <- scan("list3.txt")
# late_patho.list <- scan("list2.txt")
# patho.list <- scan("list4.txt")
# prot.combined2$patho <- recode(prot.combined2$patients, "late_patho.list=late_patho; early_patho.list=early_patho; else=no_patho")
# prot.combined2$disease <- recode(prot.combined2$patients, "patho.list=patho; else=no_patho")
# prot.combined2$cohort <- "MIT"
# 
# data.list[[22]] <- prot.combined2
# rm(prot.combined2)
# gc()
# 
# #############  AD brain 3 - Morabito
# 
# setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/snRNAseq_morabito/Seurat/")
# 
# prot.data <- Read10X_h5(filename = "../preprocessed/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
# prot.combined.meta <- read.csv("../preprocessed/GSE174367_snRNA-seq_cell_meta.csv", header = T, row.names = 1)
# 
# 
# prot.combined3 <- CreateSeuratObject(counts = prot.data, meta.data = prot.combined.meta,  min.cells = 3, min.features = 200) %>%
#   PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
#   subset(subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5) %>%
#   NormalizeData(normalization.method = "LogNormalize") %>%
#   FindVariableFeatures(selection.method = "vst") %>% ScaleData()
# 
# prot.combined3$patient <- prot.combined3$SampleID
# prot.combined3$disease <- prot.combined3$Diagnosis
# no_patho <- "C"
# patho <- "AD"
# patho.list <- c("AD")
# prot.combined3$disease <- recode(prot.combined3$disease, "patho.list=patho; else=no_patho")
# prot.combined3$cohort <- "UCI"
# 
# rm(prot.data)
# rm(prot.combined.meta)
# gc()
# 
# data.list[[23]] <- prot.combined3
# 
# rm(prot.combined3)
# gc()

################ integrated all
########  Merge everything
setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/AD1.2.5/Seurat/")

files.set <- c(
  "AD/GSM4775561_AD1/",
  "AD/GSM4775562_AD2/",
  "AD/GSM4775563_AD4/",
  "AD/GSM4775564_AD5/",
  "AD/GSM4775565_AD6/",
  "AD/GSM4775566_AD8/",
  "AD/GSM4775567_AD9/",
  "AD/GSM4775568_AD10/",
  "AD/GSM4775569_AD13/",
  "AD/GSM4775570_AD19/",
  "AD/GSM4775571_AD20/",
  "AD/GSM4775572_AD21/",
  "C/GSM4775573_NC3/",
  "C/GSM4775574_NC7/",
  "C/GSM4775575_NC11/",
  "C/GSM4775576_NC12/",
  "C/GSM4775577_NC14/",
  "C/GSM4775578_NC15/",
  "C/GSM4775579_NC16/",
  "C/GSM4775580_NC17/",
  "C/GSM4775581_NC18/",
  "prot.combined2",
  "prot.combined3"
)

# Name layers
names(data.list) <- files.set

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

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F,
        #ncol = 2,
        label = T,
        repel = T,
        group.by = "patient"
        #split.by = "cohort2"
)#, combine = F)


## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated_3000feats.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated_3000feats.Rds")


## Dimensionality reduction of integrated data
samps <- read.csv("samples.csv", header = T)
samples <- t(samps$Var1)
prot.combined <- subset(prot.combined, subset = patient %in% samples)
prot.combined <- RunHarmony(prot.combined, group.by.vars = c("patient","cohort"))
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

DimPlot(prot.combined, reduction = "umap", raster = F, 
        #ncol = 3,
        label = T,
        repel = T,
        #group.by = "seurat_clusters",
        #split.by = "cohort"
)#, combine = F)

## Save and Load data
saveRDS(object = prot.combined, file = "obj_harmony_patient.cohort_reg.out.mito.ncounts.annotated.v3.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_harmony_patient.cohort_reg.out.mito.ncounts.annotated.v3.Rds")

myeloids <- subset(prot.combined, idents = c("75","68","59","70",#"76",
                                             "19","7"))
cd8.pos <- subset(prot.combined, subset = (CD8A > 0 | CD8B > 0) & (CD3E > 0 | CD3D > 0))
cd45.pos <- subset(prot.combined, subset = PTPRC > 0)

## Save and Load data
saveRDS(object = myeloids, file = "myeloids_harmony_patient.cohort_reg.out.mito.ncounts.v3.Rds")
#rm(prot.combined)
myeloids <- readRDS("./myeloids_harmony_patient.cohort_reg.out.mito.ncounts.v3.Rds")

DotPlot(prot.combined, features = c( "GFAP","SLC1A3","AQP4","LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4",   #Astrocytes
                                     "FLT1","CLDN5", "VTN","ITM2A", "VWF", "FAM167B","BMX","CLEC1B",    #Endothelial_cells
                                     "SLC17A6","SLC17A7","NRGN","CAMK2A", "SATB2", "COL5A1","SDK2","NEFM",    #Excitatory_neurons
                                     "SLC32A1","GAD1","GAD2","TAC1","PENK","SST","NPY","MYBPC1","PVALB","GABBR2",   #Inhibitory_neurons
                                     "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",   #Microglia
                                     "OLIG2", "MBP","MOBP","PLP1","MOG","CLDN11","MYRF","GALC","ERMN","MAG",   #Oligodendrocytes
                                     "VCAN","CSPG4","PDGFRA", "SOX10","NEU4", "PCDH15","GPR37L1","C1QL1","CDO1","EPN2",   #Oligodendrocyte_precursor_cells
                                     "AMBP","HIGD1B","COX4I2", "AOC3","PDE5A","PTH1R","P2RY14","ABCC9","KCNJ8","CD248",  #Pericytes
                                     "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2", #Macrophages
                                     "CD8A", "CD3E","PTPRC"
),
#cols = c("blue","blue","blue"),#"green","yellow","gray","pink","brown","lightblue"), 
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
cluster.idents = T, #group.by = "seurat_clusters",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()

###### cell type annotation

prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Oli_1",
                              `1`= "Oli_2",
                              `2`= "Opc_1",
                              `3`= "Ast_1",
                              `4`= "Oli_3",
                              `5`= "Oli_4",
                              `6`= "Ex_1",
                              `7`= "Mic_1",
                              `8`= "Ast_2",
                              `9`= "Oli_5",
                              `10`= "Ex_2",
                              `11`= "Ex_3",
                              `12`= "In_1",
                              `13`= "Oli_6",
                              `14`= "Ex_4",
                              `15`= "Oli_7",
                              `16`= "Oli_8",
                              `17`= "In_2",
                              `18`= "Oli_9",
                              `19`= "Mic_2",
                              `20`= "Oli_10",
                              `21`= "Ex_5",
                              `22`= "In_3",
                              `23`= "Oli_11",
                              `24`= "In_4",
                              `25`= "Oli_12",
                              `26`= "Ex_6",
                              `27`= "Ex_7",
                              `28`= "In_5",
                              `29`= "Oli_13",
                              `30`= "Ex_8",
                              `31`= "Ast_3",
                              `32`= "Oli_14",
                              `33`= "In_6",          ###  SDK2+ CD36+
                              `34`= "Ast_4",
                              `35`= "Ex_9",
                              `36`= "Ex_10",
                              `37`= "Opc_2",
                              `38`= "Ex_23",              #### SDK2+
                              `39`= "Ex_11",
                              `40`= "Ast_5",
                              `41`= "In_8",       #########
                              `42`= "Ex_12",
                              `43`= "Ex_13",
                              `44`= "NE_1",       # Neuroendocrine cells - previously identified as vSMC     ### high in Ex markers / SDK2+ CD36+
                              `45`= "Endo_1",
                              `46`= "In_7",
                              `47`= "Opc_3",
                              `48`= "Ex_14",
                              `49`= "Ex_15",          ######
                              `50`= "Ex_16",
                              `51`= "Peri_1",
                              `52`= "Ast_6",
                              `53`= "Ex_17",
                              `54`= "Oli_15",
                              `55`= "Ex_18",
                              `56`= "Ex_19",
                              `57`= "Opc_4",
                              `58`= "Opc_5",
                              `59`= "Mic_3",
                              `60`= "Ast_7",
                              `61`= "Oli_16",
                              `62`= "Opc_6",
                              `63`= "Oli_17",
                              `64`= "Opc_7",
                              `65`= "Ast_9",
                              `66`= "Opc_8",
                              `67`= "Ex_20",
                              `68`= "Mic_4",
                              `69`= "Ast_8",
                              `70`= "Mic_5",
                              `71`= "Ast_10",
                              `72`= "Opc_9",
                              `73`= "Oli_18",
                              `74`= "Ex_21",
                              `75`= "Mic_6",
                              `76`= "Opc_10",
                              `77`= "Ex_22"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "seurat_clusters"


##############

DotPlot(prot.combined, features = c( "P2RY12","CSF1R","TMEM119","CX3CR1","TREM2","SALL1","CD74","C3","CST3","HEXB","C1QA","SLC2A5","AIF1",   #Microglia
                                     "CD14","FCGR3A","FCGR1A","CD68","LY6G6C","PLAC8","CEBPB","PLAUR","SELL","CCR2","ITGAM","HP","GDA","EMILIN2","TFRC","CCR5"#, #Macrophages
                                     #"IL15", "PTPRC" 
),
col.max = 20,
cols = c("blue","blue"),
dot.scale = 10, 
cluster.idents = F, #group.by = "patho",
idents = c("75","68","59","70","76","19","7"),
#scale = F,
split.by = "disease",
group.by = "seurat_clusters"
) + RotatedAxis()

prot.combined <- BuildClusterTree(prot.combined, assay = "RNA", reduction = "umap", reorder = T)

markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = TRUE)
write.xlsx(as.data.frame(markers),rowNames = T,file="all.clusters.markers.xlsx")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

prot.combined <- ScaleData(prot.combined, assay = "RNA", features = top5$gene)
# p <- DoHeatmap(prot.combined, assay = "RNA", features = top5$gene, size = 4
# ) + NoLegend() #+ theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
# p

DotPlot(prot.combined, features = unique(top5$gene),
        col.max = 20, 
        dot.scale = 10, 
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


genes1 <- c("PTPRC","IL15","IL15RA","IFNG-AS1","CD36","PSEN1","APP","HTR2C","GDNF-AS1","MAPT","BACE1","TGFB1","TGFBR2","TGFBR1")

for (i in 1:length(genes1)){
  p1 <- FeaturePlot(prot.combined, features = genes1[i], 
                    reduction = "umap",
                    raster = F
                    #,max.cutoff = 1
                    #, min.cutoff = 1.5
                    , split.by = "disease"
  )
  ggsave(filename = paste0("feat_",genes1[i],"_by.disease.pdf"),
         plot = p1,
         width = 10,
         height = 5,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}
genes1 <- c("JAK1","IL15","IL15RA","JAK2","CD36","IRF2","CREB1","STAT1","STAT2","STAT3","IRF3","TGFB1","TGFBR2","TGFBR1",
            "TNF","IFNGR1","IFNGR2","TLR2","TLR4","MAPK14","PIAS1","PPARG","SPI1")

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "disease",
                idents = "Mic"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])))
  p1$layers[[2]]$aes_params$alpha <- 0.6
  ggsave(filename = paste0("vln_",genes1[i],"_mean_mic_all_by.disease.pdf"),
         plot = p1,
         width = 6,
         height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

FeaturePlot(prot.combined, features = c("IL16")
            #)#c("IFIT3","IFNG")#platelets
            , reduction = "umap",
            raster = F
            ,max.cutoff = 2
            #, min.cutoff = 1.5
            
            , split.by = "disease"
) + theme(legend.position = "right")

#prot.combined$cellsubtype.AD <- paste(prot.combined$cellsubtype, prot.combined$ADdiag2types, sep = "_")
p1 <- VlnPlot(prot.combined, features = c("ZIC2"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
              #split.by = "disease",
              pt.size = 0.05,
              raster = F,
              #ncol = 1,
              group.by = "disease",
              #assay = "RNA",
              #slot = "data",
              #add.noise = F,
              #log = T,
              #sort = "increasing",
              idents = #c("8","25")
              "44"
              #idents = c("0","1","3","5","8","16","18","21","22","27")#"37","55","9","16","31")#"51_Izzy","51_Public-1","51_Public-2","51_Public-3")#,"53","28","26")
) #+ #geom_boxplot(width=0.1, color="black", alpha=0.2) +
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["ZIC1"]])))
p1$layers[[2]]$aes_params$alpha <- 0.6
p1

## random subset

Idents(prot.combined) <- "orig.ident"
object_subset <- subset(prot.combined, cells = Cells(prot.combined[["RNA"]]), downsample = 10000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "RNA"
Idents(object_subset) <- "celltypes"
object_subset <- BuildClusterTree(object_subset, assay = "RNA", reduction = "umap", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "RNA", only.pos = TRUE)
write.xlsx(as.data.frame(markers), rowNames = T, file="wilcox_clus_all_markers_random_subset.xlsx")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "RNA", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "RNA", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
p

DotPlot(object_subset, features = unique(top5$gene),
        col.max = 20, 
        dot.scale = 10, 
        #cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

# ggsave(filename = "prot.combined_heat_top10_HD.png",
#        plot = p,
#        width = 30,
#        height = 30,
#        #dpi = 600,
#        device = "png")

### celltypes records
library(ggstatsplot)
prot.combined$disease.patient <- paste(prot.combined$disease, prot.combined$patient, sep = "_")
num.cells <- as.data.frame(table(prot.combined$disease.patient, prot.combined$disease.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(prot.combined$disease.patient, prot.combined$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.celltypes_disease_patient_prot.combined.csv")

df_freqs <- read.csv("freq.num.celltypes_disease_patient_prot.combined.csv",header = T)
plt <- ggbetweenstats(data = df_freqs,
                      x = disease_age,
                      y = Effector.CD8..T.cells,
                      p.adjust.method = "none",
                      type = "np"
                      
)
ggsave(filename = "Effector.CD8_freq.num_cells_np_disease_age2-vlnplot.pdf",
       plot = plt,
       width = 8,
       height = 8,
       device = "pdf")

###### DEGs

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "cluster.disease"
zk.response0 <- FindMarkers(prot.combined, ident.1 = c("44_AD"), 
                            ident.2 = c("44_C"),
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
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_ADxHC_prot.combined_cluster44_DEGs.xlsx")
rm(zk.response0)


##########   MAST with regression comparing AD x C

celltypes <- c(
  "Astrocytes",
  "Immune cells (CD45+)",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons"
)

prot.combined$celltype.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "celltype.disease"

# MAST with cohort regression
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("cohort"),
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_prot.combined_lat.vars.cohort_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

# MAST without cohort regression
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              #latent.vars = c("cohort"),
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_prot.combined_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

# Wilcoxon
# by cell types
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0.0,
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
  #write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("no230_MAST_lat.vars.disease.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("wilcox_prot.combined_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

# by clusters
celltypes <- as.character(c(0:77))
Idents(prot.combined) <- "seurat_clusters"
prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "cluster.disease"

for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0.0,
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
  #write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("no230_MAST_lat.vars.disease.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("wilcox_clusters_prot.combined_ADxHC_clus",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

####   GSEA

res <- read.xlsx("wilcox_prot.combined_ADxHC_Excitatory neurons_DEGs.xlsx")
colnames(res)[1] <- "SYMBOL"

# deg_df <- res[res$p_val_adj < 0.1 & res$avg_log2FC > 0,]

# deg_df <- res[res$p_val < 0.05,]
deg_df <- res[res$p_val_adj < 0.1 & res$avg_log2FC > 0 & res$pct.1 > 0.3,]

# Prepare ranked gene list
gene_list <- deg_df$avg_log2FC
names(gene_list) <- deg_df$SYMBOL

# Sort gene list in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

# # Convert gene symbols to Entrez IDs
# entrez_ids <- bitr(names(gene_list), 
#                    fromType = "SYMBOL", 
#                    toType = "ENTREZID", 
#                    OrgDb = org.Hs.eg.db)

# Create filtered and sorted gene list with Entrez IDs
gene_list_entrez <- gene_list[entrez_ids$SYMBOL]
names(gene_list_entrez) <- entrez_ids$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# Perform GSEA using KEGG pathways
gsea_result <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",           # human (change for other species)
  keyType = "ncbi-geneid",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

# Perform GSEA using Gene Ontology
gsea_result <- gseGO(
  geneList = gene_list,
  ont = "BP",                  # Biological Process (can use "MF" or "CC")
  OrgDb = org.Hs.eg.db,        # Organism database
  keyType = "SYMBOL",          # Gene identifier type
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 0,
  verbose = FALSE
) 

# Convert symbols to Entrez IDs
entrez_ids <- bitr(names(gene_list), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
gene_list_entrez <- gene_list[entrez_ids$SYMBOL]
names(gene_list_entrez) <- entrez_ids$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# Get MSigDB gene sets (example: Hallmark)
msigdb_sets <- msigdbr(species = "Homo sapiens", category = "H")
# Alternative: Load GMT file (download from https://www.gsea-msigdb.org)
# msigdb_sets <- read.gmt("h.all.v2023.1.Hs.entrez.gmt")

# Perform GSEA with MSigDB
gsea_result <- GSEA(
  geneList = gene_list_entrez,
  TERM2GENE = msigdb_sets[, c("gs_name", "entrez_gene")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  verbose = FALSE
)

# Visualize results
# 1. Dot plot of enriched pathways
dotplot(gsea_result, showCategory=10, split=".sign") + 
  ggtitle("GSEA - Pathway Enrichment")

dotplot(gsea_result, showCategory=10, split=".sign") + 
  ggtitle("GSEA - Gene Ontology Enrichment (BP)")

# 2. GSEA plot for top pathway
gseaplot2(gsea_result, 
          geneSetID = 1, 
          title = gsea_result$Description[1],
          pvalue_table = TRUE)

# 3. Heatmap-style plot
heatplot(gsea_result, foldChange = gene_list_entrez, showCategory=5)

# Save results
write.csv(gsea_result@result, "gsea_results.csv", row.names = FALSE)


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
setwd("/media/patrick/JANELSO/Bioinfo/ref_data/nichenet_ref/human/")
lr_network = readRDS("lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/AD1.2.5/Seurat/")

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

Idents(prot.combined) <- "seurat_clusters"

# ### circus plot attempt
# ## iNK ---> Classical Monocytes
# ## regular pipe
# nichenet_output = nichenet_seuratobj_aggregate(
#   seurat_obj = prot.combined, 
#   receiver = c("50","47","29","7","68","71","43","69","44"), 
#   condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
#   sender = c("21","41","62","20","67","57"#,"65"
#   ), 
#   ligand_target_matrix = ligand_target_matrix,
#   lr_network = lr_network,
#   weighted_networks = weighted_networks,
#   top_n_ligands = 130)
# 
# ## assign cluster
# ligand_type_indication_df <- assign_ligands_to_celltype(prot.combined,
#                                                         nichenet_output$top_ligands,
#                                                         celltype_col = "celltypes") 
# 
# ligand_type_indication_df %>% head()
# ligand_type_indication_df$ligand_type %>% table()
# 
# head(nichenet_output$ligand_target_df)
# active_ligand_target_links_df <- nichenet_output$ligand_target_df
# active_ligand_target_links_df$target_type <- "AD-DE" # needed for joining tables
# circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
#                                            active_ligand_target_links_df,
#                                            cutoff = 0.40) 
# 
# head(circos_links)
# 
# 
# ligand_colors <- c("General" = "#377EB8", "Exhausted CD8+ T cells" = "#4DAF4A", 
#                    "Plasma B cells" = "#984EA3",
#                    "Platelets" = "#FF7F00"#, 
#                    #,"DC" = "#FFFF33"
#                    , "pDC" = "#F781BF"
#                    #,"CD8 T"= "#E41A1C"
# ) 
# target_colors <- c("AD-DE" = "#999999") 
# 
# vis_circos_obj <- prepare_circos_visualization(circos_links,
#                                                ligand_colors = ligand_colors,
#                                                target_colors = target_colors,
#                                                celltype_order = NULL) 
# make_circos_plot(vis_circos_obj, transparency = TRUE,  args.circos.text = list(cex = .5)) 

####### Analysis by cell types
Idents(prot.combined) <- "celltypes"

## NK ---> Monocytes
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = prot.combined, 
  receiver = c("CD8+ T cells"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("Microglia M0","Microglia MgND","Monocytes / Macrophages"
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
receiver = effect.cd8
expressed_genes_receiver = get_expressed_genes(receiver, prot.combined, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = mono
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
best_upstream_ligands = ligand_activities %>% top_n(126, aupr_corrected) %>% 
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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="all_mono_2_effect.CD8_nichenet_vis_ligand_target.xlsx")
# vis_ligand_target = vis_ligand_target[1:113,1:140]
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                                                                    color = "purple",legend_position = "top", y_axis = T,
                                                                    x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic", size = 4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.125,0.25))#c(0,0.0045,0.0090))
p_ligand_target_network ###  Plot!!!!

file1 <- read.xlsx("all_mono_2_effect.CD8_nichenet_vis_ligand_target3.xlsx")
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
  ggtitle("Monocytes -> Effector CD8+ cells") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Purples", direction = 1, limits = c(0,1)* max(abs(file1$mean.interaction)), 
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


















###################################        Myeloids analysis          ##################################
#################
 
### Normalize and scale merged obj
myeloids <- NormalizeData(myeloids, normalization.method = "LogNormalize")
myeloids <- FindVariableFeatures(myeloids, selection.method = "vst", nfeatures = 1000)
myeloids <- ScaleData(myeloids, vars.to.regress = c("cohort","patient"), #latent.data = "nFeature_RNA", 
                           model.use = "linear")#, features = all.genes)
gc()

### Dimensionality reduction and integration
myeloids <- RunPCA(myeloids, npcs = 50)
gc()
ElbowPlot(myeloids, ndims = 50)
myeloids <- FindNeighbors(myeloids, dims = 1:30, reduction = "pca")
myeloids <- FindClusters(myeloids, resolution = 0.9, cluster.name = "unintegrated_clusters"#,algorithm = "leiden"
)
myeloids <- RunUMAP(myeloids, #umap.method = "umap-learn", 
                         dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()

DimPlot(myeloids, reduction = "umap.unintegrated", raster = F,
        ncol = 2,
        label = T,
        repel = T,
        #group.by = "patient"
        split.by = "cohort"
)#, combine = F)


## Save and Load data
saveRDS(object = myeloids, file = "myeloids_harmony_patient.cohort_reg.out.mito.ncounts.annotated.v3.Rds")
# #rm(myeloids)
myeloids <- readRDS("./myeloids_harmony_patient.cohort_reg.out.mito.ncounts.annotated.v3.Rds")

## Dimensionality reduction of integrated data
myeloids <- RunHarmony(myeloids, group.by.vars = c("patient","cohort"))
gc()
ElbowPlot(myeloids, ndims = 50, reduction = "harmony")
myeloids <- RunUMAP(myeloids, dims = 1:20, reduction = "harmony", reduction.name = "umap")
myeloids <- FindNeighbors(myeloids, reduction = "harmony", dims = 1:20)
myeloids <- FindClusters(myeloids, resolution = 0.4, cluster.name = "harmony_clusters")
gc()

DimPlot(myeloids, reduction = "umap", raster = F, 
        ncol = 2,
        label = T,
        repel = T,
        #group.by = "seurat_clusters",
        split.by = "disease"
)#, combine = F)

DotPlot(myeloids, features = c(     #"TNF","IFNG",
  "KLRC1","NCAM1","IL2RB", # iNKs
  "PRF1","GZMM","GZMH","GZMA", # mNKs
  "CD3E","CD3D", # T cells
  "CD8A","CD8B","CCL5","CD244",#"PTPRC",
  "CD4","S100A4","SELL", # T cells
  "IL17F","ARG1","IL17RA",
  #"TRDC","TRGC2","TRGV9","TRGC1", #gamma delta
  "ITGB2","PECAM1","IL3RA","LAMP1", # pDC
  "NRP1","ITGAX","FCER1A","CCR7", # DC
  "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","GDA","EMILIN2", #Macrophages
  "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB","C1QA","CX3CR1","TMEM119","TREM2","SLC2A5" # microglia
  #,"CD19","IGKC","IGHM","CD27","CD1D","CD22"#,"CD86","MS4A1","IGLC2","IGLC3","IGLL5" # B cells
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()

DotPlot(myeloids, features = c(      "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",   #Microglia
                                     "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2" #Macrophages
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()



myeloids <- BuildClusterTree(myeloids, assay = "RNA", reduction = "umap", reorder = T)

markers <- FindAllMarkers(myeloids, assay = "RNA", only.pos = TRUE)
markers <- markers[markers$p_val_adj < 0.05,]
#write.xlsx(as.data.frame(markers),rowNames = T,file="myeloid.clusters.markers.xlsx")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

myeloids <- ScaleData(myeloids, assay = "RNA", features = top5$gene)
p <- DoHeatmap(myeloids, assay = "RNA", features = top5$gene, size = 4
) + NoLegend() #+ theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
p


DotPlot(myeloids, features = unique(top5$gene),
        col.max = 20, 
        dot.scale = 10, 
        cluster.idents = F, #group.by = "patho",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


myeloids$cluster.disease <- paste(myeloids$seurat_clusters, myeloids$disease, sep = "_")
Idents(myeloids) <- "cluster.disease"

DotPlot(myeloids, features = c(#"CD8B","CD3E",#"IL2RB","GSDMB","KLRB1",#"LTB","PDCD1","FASLG","IFNG","TNF","BATF",#"IL7R"
  "CD8A","IRF2","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP" 
),
col.max = 20,
#cols = c("blue","blue"),
dot.scale = 10, 
cluster.idents = F, #group.by = "patho",
idents = c("5_AD","5_C"),
scale = T,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

Idents(myeloids) <- "seurat_clusters"

clus1.markers <- c(
  #"STARD13",
  #"CD163",
  "TPRG1",
  #"ZNF804A",
  "PLXND1",
  "CD163L1",
  "PPARG",
  #"MCTP1",
  "LINC00278",
  #"IQGAP1",
  "DIRC3",
  "IFI44L",
  "AC004980.1",
  "MRTFA",
  #"NBPF1",
  #"LITAF",
  "CMAHP",
  #"MAFB",
  "ESR1",
  "TNFAIP2",
  "PLXNC1",
  "DOP1B",
  "AC009120.2",
  "PLA2G7",
  #"MS4A4E",
  #"CPNE8",
  "FAM20A",
  "GRAMD4",
  "MARCO","MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1"
)
DotPlot(myeloids, features = unique(clus1.markers),
col.max = 20,
#cols = c("blue","blue"),
dot.scale = 10, 
cluster.idents = F, #group.by = "patho",
idents = c("9_AD","9_C"),
#scale = F,
#split.by = "disease"
) + RotatedAxis() + coord_flip()


DotPlot(myeloids, features = c(#"MARCO","MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1",
  #"CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7"
  "MARCO","MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2","CD4", #"GDA",
  "CD8A","CD8B","CD3E","CD3D","IL2RB","KLRB1","IL7R","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TRAC","TRBC1"#,"TRDC","TRGC1"
                               ),
        col.max = 20,
        #cols = c("blue","blue"),
        dot.scale = 10, 
        cluster.idents = F, #group.by = "patho",
        #idents = c("9_AD","9_C"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

FeaturePlot(myeloids, features = c("IL1RAPL1"), 
            reduction = "umap",
            #max.cutoff = 1,
            raster = F,
            ncol = 2,
            split.by = "disease"
            )

FeaturePlot(myeloids, features = c("TREM2","TMEM119","CLEC7A",""), 
            reduction = "umap",
            #max.cutoff = 1,
            raster = F,
            ncol = 2)

## DEGs analysis

##########   MAST with regression comparing AD x C

celltypes <- c(
  "Microglia M0",
  "Microglia MgND",
  "Monocytes . Macrophages",
  "CD8+ T cells"
)

myeloids$celltype.disease <- paste(myeloids$celltypes2, myeloids$disease, sep = "_")
Idents(myeloids) <- "celltype.disease"

# MAST with cohort regression
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              latent.vars = c("cohort"),
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_myeloids_lat.vars.cohort_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

# MAST without cohort regression
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "MAST",
                              min.pct = 0.0,
                              min.diff.pct = -Inf,
                              verbose = TRUE,
                              only.pos = FALSE,
                              max.cells.per.ident = Inf,
                              random.seed = 1,
                              #latent.vars = c("cohort"),
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_myeloids_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

# Wilcoxon
for (i in 1:length(celltypes)) {
  zk.response0 <- FindMarkers(myeloids, ident.1 = paste0(celltypes[i],"_AD"),
                              ident.2 = paste0(celltypes[i],"_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0.0,
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
  #write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("no230_MAST_lat.vars.disease.age2.sex.cohort.Harvard_pct0.3_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("wilcox_myeloids_ADxHC_",celltypes[i],"_DEGs.xlsx"))
  rm(zk.response0)
  gc()
}

## by cluster
myeloids$cluster.disease <- paste(myeloids$seurat_clusters, myeloids$disease, sep = "_")
Idents(myeloids) <- "cluster.disease"
zk.response0 <- FindMarkers(myeloids, ident.1 = c("1_AD"), 
                            ident.2 = c("1_C"),
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
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_ADxHC_myeloids_cluster1_DEGs.xlsx")
rm(zk.response0)

####   GSEA

res <- read.xlsx("wilcox_myeloids_ADxHC_Microglia MgND_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_Microglia M0_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_CD8+ T cells_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_Monocytes . Macrophages_DEGs.xlsx")
colnames(res)[1] <- "SYMBOL"

# deg_df <- res[res$p_val_adj < 0.1 & res$avg_log2FC > 0,]

deg_df <- res[res$p_val < 0.05,]
deg_df <- res[res$p_val_adj < 0.1,]

# Prepare ranked gene list
gene_list <- deg_df$avg_log2FC
names(gene_list) <- deg_df$SYMBOL

# Sort gene list in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(names(gene_list), 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Create filtered and sorted gene list with Entrez IDs
gene_list_entrez <- gene_list[entrez_ids$SYMBOL]
names(gene_list_entrez) <- entrez_ids$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

# Perform GSEA using KEGG pathways
gsea_result <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",           # human (change for other species)
  keyType = "ncbi-geneid",
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE,
  eps = 0
)

# Visualize results
# 1. Dot plot of enriched pathways
dotplot(gsea_result, showCategory=10, split=".sign") + 
  ggtitle("GSEA - Pathway Enrichment")

# 2. GSEA plot for top pathway
gseaplot2(gsea_result, 
          geneSetID = 1, 
          title = gsea_result$Description[1],
          pvalue_table = TRUE)

# 3. Heatmap-style plot
heatplot(gsea_result, foldChange = gene_list_entrez, showCategory=5)

# Save results
write.csv(gsea_result@result, "gsea_results.csv", row.names = FALSE)








## plots

Idents(myeloids) <- "seurat_clusters"
p1 <- VlnPlot(myeloids, features = c("IL1RAPL1"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
              #split.by = "disease",
              pt.size = 0.05,
              raster = F,
              #ncol = 1,
              group.by = "disease",
              #assay = "RNA",
              #slot = "data",
              #add.noise = F,
              #log = T,
              #sort = "increasing",
              idents = c("4")
              #idents = c("0","1","3","5","8","16","18","21","22","27")#"37","55","9","16","31")#"51_Izzy","51_Public-1","51_Public-2","51_Public-3")#,"53","28","26")
) + #geom_boxplot(width=0.1, color="black", alpha=0.2) +
stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
scale_y_continuous(limits = c(0.000, 6))
p1$layers[[2]]$aes_params$alpha <- 0.6
p1

genes1 <- c("PTPRC","IL15","IL15RA","IFNG-AS1","CD36","PSEN1","APP","HTR2C","GDNF-AS1","MAPT","BACE1","TGFB1","TGFBR2","TGFBR1","APOE")
genes1 <- c("CLEC7A","APOE","SIGLEC1","LGALS3","GPX3","CCL2","IL1B","CXCL16","IRF8","TGFA")

for (i in 1:length(genes1)){
  p1 <- FeaturePlot(myeloids, features = genes1[i], 
                    reduction = "umap",
                    raster = F
                    #,max.cutoff = 1
                    #, min.cutoff = 1.5
                    #, split.by = "disease"
  ) & theme(legend.position = "right")
  ggsave(filename = paste0("feat_myeloids_",genes1[i],".pdf"),
         plot = p1,
         width = 5,
         height = 5,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}
genes1 <- c("JAK1","IL15","IL15RA","JAK2","CD36","IRF2","CREB1","STAT1","STAT2","STAT3","IRF3","TGFB1","TGFBR2","TGFBR1",
            "TNF","IFNGR1","IFNGR2","TLR2","TLR4","MAPK14","PIAS1","PPARG","SPI1")
genes1 <- c("CD8A","IRF2","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP")

for (i in 1:length(genes1)){
  p1 <- VlnPlot(myeloids, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "disease",
                idents = c("2")
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])))
  p1$layers[[2]]$aes_params$alpha <- 0.6
  ggsave(filename = paste0("vln_",genes1[i],"_mean_clus2_myeloids_by.disease.pdf"),
         plot = p1,
         width = 6,
         height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

genes1 <- c("JAK1","IL15","CD36","STAT3","IRF3","TGFB1","TGFBR1",
            "IFNGR2","TLR2","MAPK14",
            "PIAS1","PPARG")
genes1 <- c("PTPRC","NCAM1","RORA","TIGIT","GPX3","CCL2","IL1B","CXCL16","IRF8","TGFA")
genes1 <- c("ITGAX","CD83","CD1C","NRP1","CLEC4C","CD86","IL3RA","CD80","CD1A","CD40")
genes1 <- c("IL1RAPL1","KIR3DL1","LAG3","HAVCR2","KLRG1","CD96","KLRC1")
DotPlot(myeloids, features = genes1,
        col.max = 20,
        #cols = c("blue","blue"),
        dot.scale = 10, 
        cluster.idents = T, #group.by = "patho",
        #idents = c("0_AD","0_C","1_AD","1_C","2_AD","2_C"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


myeloids <- RenameIdents(myeloids, 
                              `0`= "Microglia_1",
                              `1`= "Monocytes / Macrophages_1",
                              `2`= "MgND_1",
                              `3`= "Exhausted NK_1",
                              `4`= "Exhausted NK_2",
                              `5`= "Exhausted NK_3",
                              `6`= "Microglia_2",
                              `7`= "Microglia_3",
                              `8`= "Exhausted NK_4",
                              `9`= "Effector CD8+ T cells_1",
                              `10`= "Microglia_4")


myeloids$celltypes <- sub("(.*)_.*", "\\1", Idents(myeloids))
Idents(myeloids) <- "celltypes"
Idents(myeloids) <- "seurat_clusters"

### celltypes records
library(ggstatsplot)
myeloids$disease.patient <- paste(myeloids$disease, myeloids$patient, sep = "_")
num.cells <- as.data.frame(table(myeloids$disease.patient, myeloids$disease.patient))
num.cells <- num.cells[num.cells[,3] !=0,][,2:3]
num.cells.celltype <- as.data.frame.matrix(table(myeloids$disease.patient, myeloids$celltypes))
freq.num.cells.celltype <- num.cells.celltype / num.cells[,2]*100
tfreq.num.cells.celltype <- t(freq.num.cells.celltype)
write.csv(as.data.frame(freq.num.cells.celltype), file="freq.num.celltypes_disease_patient_myeloids.csv")

df_freqs <- read.csv("freq.num.celltypes_disease_patient_myeloids.csv",header = T)
plt <- ggbetweenstats(data = df_freqs,
                      x = disease,
                      y = Monocytes...Macrophages,
                      p.adjust.method = "none",
                      type = "p"
                      
)
ggsave(filename = "monocytes_freq.num_cells_p_disease-vlnplot.pdf",
       plot = plt,
       width = 4,
       height = 8,
       device = "pdf")

######################
phagocytosis <- c(
  #"HVCN1",
  #"ATP6V0D1",
  "NCF1",
  "ATP6V1G1",
  #"ATP6V1B2",
  "CYBB",
  "RAC2",
  #"FCN1",
  "CD93",
  #"LYST",
  "NCF2",
  #"ANXA1",
  "ITGAM",
  "RAB14",
  "CLCN3",
  #"PRKCD",
  #"CORO1A",
  #"SPG11",
  #"ITGAL",
  #"P2RX7",
  #"IRF8",
  #"PECAM1",
  "TICAM2",
  "FCGR1A",
  #"CCR2",
  #"NCF4",
  "CD14",
  #"ICAM3",
  #"PAK1",
  #"ARHGAP25",
  "ITGB2",
  #"VAV1",
  #"BIN2",
  "RAC1",
  #"HCK",
  "CD36",
  #"NFKB1",
  "ABL1",
  #"SPI1",
  "JUN",
  "JUNB","FOS","FOSL2",
  "MYD88"
)
my_comparisons <- c("C","AD")
myeloids$disease <- factor(myeloids$disease, levels = my_comparisons)
DotPlot(myeloids, features = phagocytosis
        ,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F,
        idents = "Microglia",
        # scale = F,
        group.by = "disease"
        #split.by = "disease"
) + RotatedAxis() + coord_flip()



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
setwd("/media/patrick/JANELSO/Bioinfo/ref_data/nichenet_ref/human/")
lr_network = readRDS("lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/Brain/AD1.2.5/Seurat/")

lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

Idents(myeloids) <- "seurat_clusters"

# ### circus plot attempt
# ## iNK ---> Classical Monocytes
# ## regular pipe
# nichenet_output = nichenet_seuratobj_aggregate(
#   seurat_obj = myeloids, 
#   receiver = c("50","47","29","7","68","71","43","69","44"), 
#   condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
#   sender = c("21","41","62","20","67","57"#,"65"
#   ), 
#   ligand_target_matrix = ligand_target_matrix,
#   lr_network = lr_network,
#   weighted_networks = weighted_networks,
#   top_n_ligands = 130)
# 
# ## assign cluster
# ligand_type_indication_df <- assign_ligands_to_celltype(myeloids,
#                                                         nichenet_output$top_ligands,
#                                                         celltype_col = "celltypes") 
# 
# ligand_type_indication_df %>% head()
# ligand_type_indication_df$ligand_type %>% table()
# 
# head(nichenet_output$ligand_target_df)
# active_ligand_target_links_df <- nichenet_output$ligand_target_df
# active_ligand_target_links_df$target_type <- "AD-DE" # needed for joining tables
# circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
#                                            active_ligand_target_links_df,
#                                            cutoff = 0.40) 
# 
# head(circos_links)
# 
# 
# ligand_colors <- c("General" = "#377EB8", "Exhausted CD8+ T cells" = "#4DAF4A", 
#                    "Plasma B cells" = "#984EA3",
#                    "Platelets" = "#FF7F00"#, 
#                    #,"DC" = "#FFFF33"
#                    , "pDC" = "#F781BF"
#                    #,"CD8 T"= "#E41A1C"
# ) 
# target_colors <- c("AD-DE" = "#999999") 
# 
# vis_circos_obj <- prepare_circos_visualization(circos_links,
#                                                ligand_colors = ligand_colors,
#                                                target_colors = target_colors,
#                                                celltype_order = NULL) 
# make_circos_plot(vis_circos_obj, transparency = TRUE,  args.circos.text = list(cex = .5)) 

####### Analysis by cell types
Idents(myeloids) <- "celltypes"

## NK ---> Monocytes
## regular pipe celltypes
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = myeloids, 
  receiver = c("CD8+ T cells"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("Microglia M0","Microglia MgND","Monocytes / Macrophages"
  ), 
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks)


###### Analysis by clusters
## NK ---> Monocytes
## regular pipe clusters
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = myeloids2, 
  receiver = c("80","18","48","19","58","81","78","50","47","29","7","68","71","43","69","44"), 
  condition_colname = "disease", condition_oi = "AD", condition_reference = "C", 
  sender = c("67","57","28","21","41","62","1","2","20"#,"65"
  ), 
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
# DotPlot(myeloids %>% subset(idents = c("10","23","38","30")), cols = c("royalblue","royalblue"),
#         features = nichenet_output$top_targets %>% rev(), split.by = "disease") + RotatedAxis()
nichenet_output$ligand_activity_target_heatmap
nichenet_output$ligand_receptor_heatmap


### cell types by clusters
eff.nkt.cd8 <- c("9","14","3","30","52","72","12","28","17","25","64")
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
receiver = effect.cd8
expressed_genes_receiver = get_expressed_genes(receiver, myeloids, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = mono
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, myeloids, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
seurat_obj_receiver= subset(myeloids, idents = receiver)
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
best_upstream_ligands = ligand_activities %>% top_n(126, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
best_upstream_ligands = file1$ligands
DotPlot(myeloids, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()  ### Plot!!!!


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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="all_mono_2_effect.CD8_nichenet_vis_ligand_target.xlsx")
# vis_ligand_target = vis_ligand_target[1:113,1:140]
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                                                                    color = "purple",legend_position = "top", y_axis = T,
                                                                    x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic", size = 4)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.125,0.25))#c(0,0.0045,0.0090))
p_ligand_target_network ###  Plot!!!!

file1 <- read.xlsx("all_mono_2_effect.CD8_nichenet_vis_ligand_target3.xlsx")
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
  ggtitle("Monocytes -> Effector CD8+ cells") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Purples", direction = 1, limits = c(0,1)* max(abs(file1$mean.interaction)), 
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
DE_table_all = Idents(myeloids) %>% levels() %>% intersect(sender_celltypes) %>% 
  lapply(get_lfc_celltype, seurat_obj = myeloids, condition_colname = "disease", 
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
rotated_dotplot = DotPlot(myeloids %>% subset(seurat_clusters %in% sender_celltypes), # subset(seurat_clusters %in% sender_celltypes), 
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



########################################################################





















# #################################  ###  #################    CD45+ cells
# 
# ### Normalize and scale merged obj
# cd45.pos <- NormalizeData(cd45.pos, normalization.method = "LogNormalize")
# cd45.pos <- FindVariableFeatures(cd45.pos, selection.method = "vst", nfeatures = 1000)
# cd45.pos <- ScaleData(cd45.pos, vars.to.regress = c("cohort","patient"), #latent.data = "nFeature_RNA", 
#                       model.use = "linear")#, features = all.genes)
# gc()
# 
# ### Dimensionality reduction and integration
# cd45.pos <- RunPCA(cd45.pos, npcs = 50)
# gc()
# ElbowPlot(cd45.pos, ndims = 50)
# cd45.pos <- FindNeighbors(cd45.pos, dims = 1:30, reduction = "pca")
# cd45.pos <- FindClusters(cd45.pos, resolution = 0.9, cluster.name = "unintegrated_clusters")
# cd45.pos <- RunUMAP(cd45.pos, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# gc()
# 
# DimPlot(cd45.pos, reduction = "umap.unintegrated", raster = F,
#         ncol = 2,
#         label = T,
#         repel = T,
#         #group.by = "patient"
#         split.by = "cohort"
# )#, combine = F)
# 
# DotPlot(cd45.pos, features = c(      "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",   #Microglia
#                                      "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2" #Macrophages
# ),
# col.max = 20, 
# dot.scale = 10, 
# cluster.idents = T, #group.by = "patho",
# #scale = F,
# #split.by = "cohort"
# ) + RotatedAxis()
# 
# ## Dimensionality reduction of integrated data
# cd45.pos <- RunHarmony(cd45.pos, group.by.vars = c("patient","cohort"))
# gc()
# ElbowPlot(cd45.pos, ndims = 50, reduction = "harmony")
# cd45.pos <- RunUMAP(cd45.pos, dims = 1:21, reduction = "harmony", reduction.name = "umap")
# cd45.pos <- FindNeighbors(cd45.pos, reduction = "harmony", dims = 1:21)
# cd45.pos <- FindClusters(cd45.pos, resolution = 0.45, cluster.name = "harmony_clusters")
# gc()
# 
# DimPlot(cd45.pos, reduction = "umap", raster = F, 
#         #ncol = 2,
#         label = T,
#         repel = T,
#         #group.by = "seurat_clusters",
#         #split.by = "cohort"
# )#, combine = F)
# 
# 
# ## Save and Load data
# saveRDS(object = cd45.pos, file = "cd45.pos_harmony_patient.cohort_reg.out.mito.ncounts.v3.Rds")
# # #rm(cd45.pos)
# # cd45.pos <- readRDS("./obj_unintegrated_3000feats.Rds")
# 
# DotPlot(cd45.pos, features = c("MARCO","MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2","CD4", #"GDA",
#                                "CD8A","CD8B","CD3E","CD3D","IL2RB","KLRB1","IL7R","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TRAC","TRBC1"#,"TRDC","TRGC1"
# ),
# col.max = 20,
# #cols = c("blue","blue"),
# dot.scale = 10, 
# cluster.idents = F, #group.by = "patho",
# #idents = c("9_AD","9_C"),
# #scale = F,
# #split.by = "disease"
# ) + RotatedAxis() + coord_flip()
# 
# FeaturePlot(cd45.pos, features = c("IL15"), 
#             #reduction = "umap.unintegrated",
#             reduction = "umap",
#             #max.cutoff = 1,
#             raster = F,
#             #ncol = 2
#             )
# 
# 
# 
# DotPlot(cd45.pos, features = c(      "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119", "TREM2","SLC2A5",   #Microglia
#                                      "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","EMILIN2", #Macrophages "GDA",
#                                      "CD8A","CD3E","IL2RB","GZMM","PRF1","KLRB1","GZMH","GZMA","GZMB","IL7R","NKG7"
# ),
# col.max = 20, 
# dot.scale = 10, 
# cluster.idents = T, #group.by = "patho",
# #scale = F,
# #split.by = "cohort"
# ) + RotatedAxis()
# 
# DotPlot(cd45.pos, features = c(     #"TNF","IFNG",
#                                     "KLRC1","NCAM1","IL2RB", # iNKs
#                                     "PRF1","GZMM","GZMH","GZMA", # mNKs
#                                     "CD3E","CD3D", # T cells
#                                     "CD8A","CD8B","CCL5","CD244",#"PTPRC",
#                                     "CD4","S100A4","SELL", # T cells
#                                     #"TRDC","TRGC2","TRGV9","TRGC1", #gamma delta
#                                     "ITGB2","PECAM1","IL3RA","LAMP1", # pDC
#                                     "NRP1","ITGAX","FCER1A","CCR7", # DC
#                                     "CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","GDA","EMILIN2", #Macrophages
#                                     "AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB","C1QA","CX3CR1","TMEM119","TREM2","SLC2A5" # microglia
#                                     #,"CD19","IGKC","IGHM","CD27","CD1D","CD22"#,"CD86","MS4A1","IGLC2","IGLC3","IGLL5" # B cells
# ),
# col.max = 20, 
# dot.scale = 10, 
# cluster.idents = T, #group.by = "patho",
# #scale = F,
# #split.by = "cohort"
# ) + RotatedAxis()
# 
# 
# cd45.pos <- BuildClusterTree(cd45.pos, assay = "RNA", reduction = "umap", reorder = T)
# 
# markers <- FindAllMarkers(cd45.pos, assay = "RNA", only.pos = TRUE, min.pct = 0.1)
# markers <- markers[markers$p_val_adj < 0.1,]
# write.xlsx(as.data.frame(markers),rowNames = T,file="cd45.pos.clusters.markers.xlsx")
# markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 5) %>%
#   ungroup() -> top5
# 
# cd45.pos <- ScaleData(cd45.pos, assay = "RNA", features = top5$gene)
# p <- DoHeatmap(cd45.pos, assay = "RNA", features = top5$gene, size = 4
# ) + NoLegend() #+ theme(axis.text = element_text(size = 5.5)) #+ NoLegend()
# p
# 
# 
# DotPlot(cd45.pos, features = unique(top5$gene),
#         col.max = 20, 
#         dot.scale = 10, 
#         cluster.idents = F, #group.by = "patho",
#         #scale = F,
#         #split.by = "disease"
# ) + RotatedAxis() + coord_flip()
# 
# 
# 
# cd45.pos$cluster.disease <- paste(cd45.pos$seurat_clusters, cd45.pos$disease, sep = "_")
# Idents(cd45.pos) <- "cluster.disease"
# 
# DotPlot(cd45.pos, features = c(#"CD8B","CD3E",#"IL2RB","GSDMB","KLRB1",#"LTB","PDCD1","FASLG","IFNG","TNF","BATF",#"IL7R"
#   "CD8A","IRF2","PRF1","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP" 
# ),
# col.max = 20,
# #cols = c("blue","blue"),
# dot.scale = 10, 
# cluster.idents = F, #group.by = "patho",
# idents = c("7_AD","7_C"),
# scale = F,
# #split.by = "disease"
# ) + RotatedAxis() + coord_flip()
# 
# Idents(cd45.pos) <- "seurat_clusters"
# 
# DotPlot(cd45.pos, features = c("CD8A","CD3E","IL2RB","GZMM","PRF1","KLRB1","GZMH","GZMA","GZMB","IL7R","NKG7"),
#         col.max = 20,
#         #cols = c("blue","blue"),
#         dot.scale = 10, 
#         cluster.idents = F, #group.by = "patho",
#         #scale = F,
#         #split.by = "disease"
# ) + RotatedAxis() + coord_flip()
# 
# FeaturePlot(cd45.pos, features = c("CD8A","CD8B","CD3D","CD3E"), 
#             #reduction = "umap.unintegrated",
#             reduction = "umap",
#             max.cutoff = 1,
#             raster = F,
#             ncol = 2)
# 
# p1 <- VlnPlot(cd45.pos, features = c("GZMB"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
#               #split.by = "disease",
#               pt.size = 0.05,
#               raster = F,
#               #ncol = 1,
#               group.by = "disease",
#               #assay = "RNA",
#               #slot = "data",
#               #add.noise = F,
#               #log = T,
#               #sort = "increasing",
#               idents = c("7")
#               #idents = c("0","1","3","5","8","16","18","21","22","27")#"37","55","9","16","31")#"51_Izzy","51_Public-1","51_Public-2","51_Public-3")#,"53","28","26")
# ) + #geom_boxplot(width=0.1, color="black", alpha=0.2) +
# stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
# scale_y_continuous(limits = c(0.000, 4))
# p1$layers[[2]]$aes_params$alpha <- 0.6
# p1
# 
# 
# cd45.pos$cluster.disease <- paste(cd45.pos$seurat_clusters, cd45.pos$disease, sep = "_")
# Idents(cd45.pos) <- "cluster.disease"
# zk.response0 <- FindMarkers(cd45.pos, ident.1 = "8_AD", 
#                             ident.2 = "8_C",
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
# #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
# write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_ADxHC_cd45.pos_cluster8_DEGs.xlsx")
# rm(zk.response0)
# 
# 
# ################### ###  ###  #################    CD8 cells
# Idents(cd8.pos) <- "disease"
# 
# DotPlot(cd8.pos, features = c("CD8B","IRF1","GZMA","PRF1","LTB","PDCD1","GZMK","GZMB","NKG7",#"IFNG","TNF",
#                               "TOX"#
#   #"GZMA","NKG7","GZMB","GZMK","PRF1","TOX"
# ),
# col.max = 20, 
# dot.scale = 10, 
# #cluster.idents = T, #group.by = "patho",
# #scale = F,
# #split.by = "cohort"
# group.by = "disease"
# ) + RotatedAxis()
# 
# zk.response0 <- FindMarkers(cd8.pos, ident.1 = "AD", 
#                             ident.2 = "C",
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
# #zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
# write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_ADxHC_cd8.pos_DEGs.xlsx")
# rm(zk.response0)
# 
# p1 <- VlnPlot(cd8.pos, features = c("GZMA"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
#               #split.by = "disease",
#               pt.size = 0.05,
#               raster = F,
#               #ncol = 1,
#               group.by = "disease",
#               #assay = "RNA",
#               #slot = "data",
#               #add.noise = F,
#               #log = T,
#               #sort = "increasing",
#               #idents = c("8","25")
#               #idents = c("0","1","3","5","8","16","18","21","22","27")#"37","55","9","16","31")#"51_Izzy","51_Public-1","51_Public-2","51_Public-3")#,"53","28","26")
# ) #+ #geom_boxplot(width=0.1, color="black", alpha=0.2) +
#   #stat_summary(fun = median, geom='point', size = 35, colour = "black", shape = 95) +
#   #scale_y_continuous(limits = c(0.000, 10))
# p1$layers[[2]]$aes_params$alpha <- 0.6
# p1
# 
# 








