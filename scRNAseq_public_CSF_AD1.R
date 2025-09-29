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
#library(bluster)
#library(pheatmap)

setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/public_data/AD/CSF/scRNAseq_CSF_AD1/Seurat/")

## Loop through h5 files and output BPCells matrices on-disk
file.dir <- "../processed/"
meta.data <- read.delim("meta.txt", header = T)
disease.state <- meta.data$condition
files.set <- meta.data$GEM

data.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i])
  data <- Read10X(data.dir = path)
  dataset_name <- files.set[i]
  mat <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 200) %>% 
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20) %>%
    SCTransform(vst.flavor = "v2", method = "glmGamPoi", vars.to.regress = "percent.mt", return.only.var.genes = F)
  mat$disease <- disease.state[i]
  mat$patient <- dataset_name
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
ElbowPlot(prot.combined, ndims = 50)
prot.combined <- FindNeighbors(prot.combined, dims = 1:30, reduction = "pca")
prot.combined <- FindClusters(prot.combined, resolution = 0.8, cluster.name = "unintegrated_clusters")
prot.combined <- RunUMAP(prot.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
gc()
DimPlot(prot.combined, reduction = "umap.unintegrated", label = T, raster=F, repel = T,
        #group.by = "condition"
        #split.by = "disease"#, ncol = 4
)
gc()

FeaturePlot(prot.combined, features = c("TOX","GZMA")
            #)#c("IFIT3","IFNG")#platelets
            , reduction = "umap.unintegrated",
            raster = F
            #,max.cutoff = 2.5
            #, min.cutoff = 1.5
            
            , split.by = "disease"
)

p1 <- VlnPlot(prot.combined, features = c("GZMA"),#c("JUN","STAT1", "CCL3", "CCL3L1"),#c("IRF1","IFNG","IFNGR1","IFNGR2"),#c("CD8A","CD4","CD19"),#c("TMEM176A","TMEM176B"), 
        #split.by = "disease",
        pt.size = 0.05,
        raster = F,
        ncol = 1,
        group.by = "disease",
        #slot = "counts",
        #add.noise = F,
        #log = T,
        #sort = "increasing",
        idents = "2"#c("2","4","16")#"37","55","9","16","31")#"51_Izzy","51_Public-1","51_Public-2","51_Public-3")#,"53","28","26")
) + #geom_boxplot(width=0.1, color="black", alpha=0.2) +
  stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000,5))
p1$layers[[2]]$aes_params$alpha <- 0.5
p1


DotPlot(prot.combined, features = c("TNF","IFNG","KLRC1","NCAM1","IL2RB", # iNKs
                                    "PRF1","GZMM","GZMH","GZMA", # mNKs
                                    "CD3E","CD3D", # T cells
                                    "CD8A","CD8B","PTPRC","CCL5",#"CD244",
                                    "CD4","S100A4","SELL", # T cells
                                    "FOXP3", "IL2RA", # Tregs
                                    "TRDC","TRDV1","TRDV2","TRGC2","TRGV9","TRGC1",#"TRGV7", #gamma delta
                                    "ITGB2","PECAM1","IL3RA","LAMP1", #pDC
                                    "CD1C","ITGAX","FCER1A","CCR7","NRP1", # DC
                                    "CD14","FCGR3A", #monocytes
                                    #"PF4", # platelets
                                    "CD19","IGKC","IGHM","CD27","CD1D","CD22","CD86","MS4A1","IGLC2","IGLC3" #, "IGLL5" b cells
),
cols = c("blue","blue","blue"),#"green","yellow","gray","pink","brown","lightblue"), 
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
#cluster.idents = T, #group.by = "patho",
#scale = F,
split.by = "disease"
) + RotatedAxis()

DotPlot(prot.combined, features = c(#"CD8A","CD8B",
                                    "IRF3","IRF2","IRF1","GZMA","GZMK","PRF1","LTB",#"PDCD1","GZMH","GZMB",
                                    "FASLG",#"IFNG","TNF",
                                    "TOX",#"BATF","TYROBP","GSDMB"
                                    "PRDX1","FCGR3A","LAG3","CD2",#"KLRF2","KIR3DL1","STAT5B",
                                    "VAMP7","KLRK1","MICA", #"CORO1A","PLEKHM2","IL18","TUBB4B","LYST","SH2D1A","CLEC2A","RAET1E","RAET1G",,"ULBP2","ULBP3","TUBB","HCST","RNF19B"
                                    "CEBPG","GRB2","UNC13D","VAV1","ARL8B","PIK3R1","PTPN6","RAB27A","SLAMF7","KIF5B"
                                    ),#cd8act,
        #cols = c("blue","blue","blue"),
        #cols = c("blue","red"),#"green","yellow","gray","pink","brown","lightblue"), 
        col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),#"Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"),
          #c("Classical Mono_1","Classical Mono_2","Intermediate Mono","Nonclassical Mono", 
          #"pDC_AD","pDC_C", 
          #"mo-DC_AD","mo-DC_C"
          #c("34","50","9","16","31","30","42"),
      idents = "NK",
        # c("NK_4_C","NK_4_AD","NK_8_C","NK_8_AD","NK_21_C","NK_21_AD",
        #  "CD8+ NKT-like_C", "CD8+ NKT-like_AD"#, "NK_AD", "NK_C"
        # c("NK_AD","NK_C","Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C","Intermediate Mono_AD","Nonclassical Mono_AD","Intermediate Mono_C","Nonclassical Mono_C"
        # ),
        dot.scale = 10,
        
        cluster.idents = T, 
        group.by = "disease",
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


## Save and Load data
saveRDS(object = prot.combined, file = "obj_unintegrated.Rds")
#rm(prot.combined)
prot.combined <- readRDS("./obj_unintegrated.Rds")



# ## integration
# prot.combined <- RunHarmony(prot.combined, group.by.vars = "sample")
# ElbowPlot(prot.combined, ndims = 50, reduction = "harmony")
# prot.combined <- RunUMAP(prot.combined, dims = 1:30, reduction = "harmony", reduction.name = "umap")
# prot.combined <- FindNeighbors(prot.combined, dims = 1:30, reduction = "harmony")
# prot.combined <- FindClusters(prot.combined, graph.name = "SCT_snn", resolution = 2.5, cluster.name = "harmony_clusters2")
# gc()
# DimPlot(prot.combined, reduction = "umap", label = T, raster=F,
#         #group.by = "condition"
#         #split.by = "condition"#, ncol = 4
# )

DotPlot(prot.combined, features = c("TNF","IFNG","TRBV20-1","TRBV5-1","CCR6","KLRC1","NCAM1","IL2RB","IL7R","TBX21","EOMES",# iNKs
                                    "PRF1","GZMM","GZMK","GZMH","GZMA","GZMB","NKG7", # mNKs
                                    "LILRB1", "KLRB1", "ZBTB16", # NKT-like
                                    "CD3E","CD3D", # T cells
                                    "CD8A","CD8B","PTPRC","CCL5", #"CD244", # T cells
                                    "CD4","S100A4","SELL", # T cells
                                    "FOXP3", "IL2RA", # Tregs
                                    "TRDC","TRDV1","TRDV2","TRGC2","TRGV9","TRGC1", #gamma delta
                                    "ITGB2","PECAM1","IL3RA","LAMP1", #pDC
                                    #"CD64", "CD68", "CD71", "CCR5","ITGAM" # Macrophages
                                    "CD1C","ITGAX","FCER1A","CCR7","NRP1", # DC
                                    "CD14","FCGR3A", #monocytes
                                    #"PF4", # platelets
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
cluster.idents = T, #group.by = "patho",
#scale = F,
#split.by = "disease"
) + RotatedAxis()

prot.markers <- FindAllMarkers(prot.combined, assay = "RNA", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.33) %>% group_by(cluster)
write.csv(as.data.frame(prot.markers), file="All_markers.csv")

prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Effector CD4+ T cells_1",
                              `1`= "Naive CD4+ T cells_1",
                              `2`= "Effector CD8+ T cells_1",
                              `3`= "NKT-like CD4+_1",
                              `4`= "Naive CD8+ T cells_1",
                              `5`= "Effector CD4+ T cells_1",
                              `6`= "NKT-like CD4+_1",
                              `7`= "Effector CD4+ T cells_1",
                              `8`= "Intermediate monocytes_1",
                              `9`= "mo-DC_1",
                              `10`= "Classical monocytes_1",
                              `11`= "NK_1",
                              `12`= "Naive CD4+ T cells_2",
                              `13`= "Memory T CD4+_1",
                              `14`= "Effector CD4+ T cells_1",
                              `15`= "pDC_1",
                              `16`= "NKT-like CD8+_1",
                              `17`= "Intermediate monocytes_2",
                              `18`= "B cells_1"
                              )


prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
prot.combined$celltypes.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")

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


genes1 <- c("PTPRC","IL15","IL15RA","IFNG-AS1","CD36","PSEN1","APP","HTR2C","GDNF-AS1","MAPT","BACE1","TGFB1","TGFBR2","TGFBR1","APOE")
genes1 <- c("CLEC7A","APOE","SIGLEC1","LGALS3","GPX3","CCL2","IL1B","CXCL16","IRF8","TGFA")

for (i in 1:length(genes1)){
  p1 <- FeaturePlot(prot.combined, features = genes1[i], 
                    reduction = "umap",
                    raster = F
                    #,max.cutoff = 1
                    #, min.cutoff = 1.5
                    #, split.by = "disease"
  ) & theme(legend.position = "right")
  ggsave(filename = paste0("feat_prot.combined_",genes1[i],".pdf"),
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
genes1 <- c("TNF","IFNG","TGFB1")

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05,
                raster = F,
                group.by = "disease",
                #idents = c("Effector CD8+ T cells")
  ) 
  p1 <- p1 + stat_summary(fun = median, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.0001, max(p1[[1]][["data"]][[genes1[i]]])))
  p1$layers[[2]]$aes_params$alpha <- 0.6
  ggsave(filename = paste0("vln_",genes1[i],"_median_global_no.0s_by.disease.pdf"),
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
DotPlot(prot.combined, features = genes1,
        col.max = 20,
        #cols = c("blue","blue"),
        dot.scale = 10, 
        cluster.idents = T, #group.by = "patho",
        #idents = c("0_AD","0_C","1_AD","1_C","2_AD","2_C"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

prot.combined <- BuildClusterTree(prot.combined, assay = "SCT", reduction = "umap.unintegrated", reorder = T)

prot.combined <- PrepSCTFindMarkers(prot.combined)
markers <- FindAllMarkers(prot.combined, assay = "SCT", only.pos = TRUE)
write.xlsx(as.data.frame(markers),rowNames = T,file="all.celltypes.markers.xlsx")
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top5

# prot.combined <- ScaleData(prot.combined, assay = "RNA", features = top5$gene)
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



######################
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
my_comparisons <- c("HC","MCI","AD")
prot.combined$disease <- factor(prot.combined$disease, levels = my_comparisons)
DotPlot(prot.combined, features = phagocytosis
        ,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F,
        idents = "Classical monocytes",
        #scale = T,
        group.by = "disease"
        #split.by = "disease"
) + RotatedAxis() + coord_flip()



############################# Identify cell types with scType

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
DimPlot(prot.combined, reduction = "umap.unintegrated",label = TRUE, repel = TRUE, 
        group.by = 'customclassif', 
        #split.by = "condition"
)   

########################################################################

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
prot.combined$celltypes.condition <- paste(prot.combined$celltypes, prot.combined$condition, sep = "_")
prot.combined$celltypes.condition <- paste(prot.combined$customclassif, prot.combined$condition, sep = "_")
prot.combined$cluster.condition <- paste(Idents(prot.combined), prot.combined$disease, sep = "_")
prot.combined$cluster.condition <- paste(prot.combined$seurat_clusters, prot.combined$condition, sep = "_")
prot.combined$celltypes.disease <- paste(prot.combined$celltypes, prot.combined$disease, sep = "_")
Idents(prot.combined) <- "celltypes.disease"
prot.combined2$celltypes <- Idents(prot.combined2)
Idents(prot.combined) <- "seurat_clusters"
Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "celltypes.condition"
Idents(prot.combined) <- "cluster.condition"
Idents(prot.combined) <- "customclassif"

clusters <- as.character(0:39)
cell.types <- c(
  "Effector CD4+ T cells",
  "Naive CD4+ T cells",
  "Effector CD8+ T cells",
  "NKT-like CD4+",
  "Naive CD8+ T cells",
  "Intermediate monocytes",
  "mo-DC",
  "Classical monocytes",
  "NK",
  "Memory T CD4+",
  "pDC",
  "NKT-like CD8+",
  "B cells"
)

for (i in 1:length(clusters)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(clusters[i], "_AD"), 
                              ident.2 = paste0(clusters[i], "_HC"),
                              test.use = "MAST", 
                              min.pct = 0.0, 
                              logfc.threshold = 0,
                              assay = "SCT", 
                              slot = "data", 
                              #latent.vars = "condition",
                              latent.vars = "nCount_RNA", 
                              pseudocount.use = 0.001,
                              verbose = FALSE
  )
  write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("MAST_SCT_ADxHC_", clusters[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(cell.types[i], "_AD"), 
                              ident.2 = paste0(cell.types[i], "_HC"),
                              slot = "data",
                              assay = "SCT",
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_SCT_ADxHC_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

# zk.response0 <- FindMarkers(prot.combined, ident.1 = "11_AD", 
#                             ident.2 = "11_HC",
#                             test.use = "MAST", 
#                             min.pct = 0.0, 
#                             logfc.threshold = 0,
#                             assay = "RNA", 
#                             slot = "counts", 
#                             latent.vars = "nCount_RNA", 
#                             #pseudocount.use = 0.001,
#                             verbose = FALSE
# )
# write.csv(as.data.frame(zk.response0), file="MAST_ADxC_clus11_DEGs.csv")
# rm(zk.response0)
# 
# 
# zk.response0 <- FindMarkers(prot.combined, ident.1 = c("2_AD","4_AD","16_AD"), 
#                             ident.2 = c("2_HC","4_HC","16_HC"),
#                             slot = "data",
#                             assay = "RNA",
#                             features = NULL,
#                             logfc.threshold = 0,
#                             test.use = "wilcox",
#                             min.pct = 0.0,
#                             min.diff.pct = -Inf,
#                             verbose = TRUE,
#                             only.pos = F,
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
# write.xlsx(as.data.frame(zk.response0), rowNames = T, file="wilcox_2.4.16_ADxHC_DEGs.xlsx")
# rm(zk.response0)
# 
# nums <- as.character(0:18)
# for (i in 1:length(nums)) {
#   zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(nums[i],"_AD"), 
#                               ident.2 = paste0(nums[i],"_HC"),
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
#   write.xlsx(as.data.frame(zk.response0), rowNames = T, file=paste0("wilcox_ADxHC_clus",nums[i],"_DEGs.xlsx"))
#   rm(zk.response0)
# }





