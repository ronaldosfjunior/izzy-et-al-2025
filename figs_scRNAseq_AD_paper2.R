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
options(Seurat.object.assay.version = "v3")
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

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/public_data/AD/Brain/AD1.2.5/Seurat/")


prot.combined <- readRDS("./obj_harmony_patient.cohort_reg.out.mito.ncounts.annotated.v3.Rds")

myeloids <- readRDS("./myeloids_harmony_patient.cohort_reg.out.mito.ncounts.v3.Rds")


################################   PPG Oleg  ###################

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/PPG_Oleg_AD/")

hav <- t(unique(read.delim("AD-GWAS.txt", header = F)))
strips <- "^\\s+|\\s+$"
hav <- gsub(strips,"",hav)
res <- read.xlsx("wilcox_prot.combined_ADxHC_Immune cells (CD45+)_DEGs.xlsx")
colnames(res)[1] <- "genes"

EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.00000215775043879316,
                FCcutoff = 0.5,
                xlim = c(-6, 6),
                ylim = c(0, 30),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 3.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Microglia - AD x C",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                maxoverlapsConnectors = Inf,#NULL
                #min.segment.length = 0.0000001,
                directionConnectors = "both",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)

ad.gwas <- res[res$genes %in% hav,]
write.xlsx(as.data.frame(ad.gwas), rowNames = T,file=paste0("microglia_AD-GWAS_genes_ADxHC.xlsx"))

myeloids$celltypes.disease <- paste(myeloids$celltypes, myeloids$disease, sep = "_")
Idents(myeloids) <- "celltypes.disease"

zk.response0 <- FindMarkers(myeloids, ident.1 = c("Microglia MgND_AD","Microglia M0_AD","Monocytes . Macrophages_AD"), 
                              ident.2 = c("Microglia MgND_C","Microglia M0_C","Monocytes . Macrophages_C"),
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
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_myeloids_ADxHC_microglia_DEGs.xlsx")
zk.response0 <- read.xlsx("wilcox_myeloids_ADxHC_microglia_DEGs.xlsx")
colnames(zk.response0)[1] <- "genes"
ad.gwas2 <- zk.response0[zk.response0$genes %in% hav,]
write.xlsx(as.data.frame(ad.gwas2), rowNames = T,file=paste0("microglia_AD-GWAS_genes_ADxHC_2.xlsx"))

Idents(myeloids) <- "celltypes"
myeloids2 <- subset(myeloids, idents = c("Microglia M0","Monocytes . Macrophages","Microglia MgND"))

avg.exp <- AverageExpression(
  myeloids2,
  assays = "RNA",
  features = hav,
  return.seurat = FALSE,
  group.by = "disease",
  add.ident = NULL,
  layer = "data",
  # slot = deprecated(),
  verbose = TRUE
)

write.xlsx(as.data.frame(avg.exp$RNA), rowNames = T,file=paste0("microglia_AD-GWAS_genes_avg.exp_2.xlsx"))

Idents(prot.combined) <- "celltypes"
avg.exp2 <- AverageExpression(
  prot.combined,
  assays = "RNA",
  features = hav,
  return.seurat = FALSE,
  group.by = "celltypes",
  add.ident = NULL,
  layer = "data",
  # slot = deprecated(),
  verbose = TRUE
)
write.xlsx(as.data.frame(avg.exp2$RNA), rowNames = T,file=paste0("brain_celltypes_all_genes_avg.exp.xlsx"))
write.xlsx(as.data.frame(avg.exp2$RNA), rowNames = T,file=paste0("brain_celltypes_AD-GWAS_genes_avg.exp.xlsx"))
hav2 <- t(read.delim("AD-GWAS_microglia.txt", header = F))

## Targets heatmap
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
heat <- pheatmap(avg.exp2$RNA[row.names(avg.exp2$RNA) %in% hav2,], cluster_rows = T, show_rownames = T, 
                 show_colnames = T, #annotation_col = sampleTable2,
                 cluster_cols = F, annotation_names_col = T, annotation_legend = T,
                 color = colors,
                 main = "AD-GWAS genes", fontsize_row = 7,
                 #labels_row = as.expression(newnames),
                 angle_col = 45,
                 scale = "row"
) 

ad.gwas3 <- zk.response0[zk.response0$genes %in% hav2,]
ad.gwas4 <- ad.gwas3[ad.gwas3$p_val < 0.05,]
ad.gwas5 <- ad.gwas4[abs(ad.gwas4$avg_log2FC) > 0.05,]
write.xlsx(as.data.frame(ad.gwas5), rowNames = T,file=paste0("microglia_AD-GWAS_ADxC_DEGs_stats.xlsx"))

EnhancedVolcano(zk.response0,
                lab = zk.response0$genes, #NA,
                selectLab = ad.gwas5$genes,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.05,
                xlim = c(-6, 6),
                ylim = c(0, 60),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Microglia - AD x C",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                maxoverlapsConnectors = Inf,#NULL
                #min.segment.length = 0.0000001,
                directionConnectors = "both",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)

###############################        fig 5     ##############################

######## fig 5 A

Idents(prot.combined) <- "seurat_clusters"
prot.combined <- RenameIdents(prot.combined, 
                              `0`= "Oligodendrocytes_1",
                              `1`= "Oligodendrocytes_2",
                              `2`= "Oligodendrocytes precursor cells_1",
                              `3`= "Astrocytes_1",
                              `4`= "Oligodendrocytes_3",
                              `5`= "Oligodendrocytes_4",
                              `6`= "Excitatory neurons_1",
                              `7`= "Immune cells (CD45+)_1",
                              `8`= "Astrocytes_2",
                              `9`= "Oligodendrocytes_5",
                              `10`= "Excitatory neurons_2",
                              `11`= "Excitatory neurons_3",
                              `12`= "Inhibitory neurons_1",
                              `13`= "Oligodendrocytes_6",
                              `14`= "Excitatory neurons_4",
                              `15`= "Oligodendrocytes_7",
                              `16`= "Oligodendrocytes_8",
                              `17`= "Inhibitory neurons_2",
                              `18`= "Oligodendrocytes_9",
                              `19`= "Immune cells (CD45+)_2",
                              `20`= "Oligodendrocytes_10",
                              `21`= "Excitatory neurons_5",
                              `22`= "Inhibitory neurons_3",
                              `23`= "Oligodendrocytes_11",
                              `24`= "Inhibitory neurons_4",
                              `25`= "Oligodendrocytes_12",
                              `26`= "Excitatory neurons_6",
                              `27`= "Excitatory neurons_7",
                              `28`= "Inhibitory neurons_5",
                              `29`= "Oligodendrocytes_13",
                              `30`= "Excitatory neurons_8",
                              `31`= "Astrocytes_3",
                              `32`= "Oligodendrocytes_14",
                              `33`= "Inhibitory neurons_6",          ###  SDK2+ CD36+
                              `34`= "Astrocytes_4",
                              `35`= "Excitatory neurons_9",
                              `36`= "Excitatory neurons_10",
                              `37`= "Oligodendrocytes precursor cells_2",
                              `38`= "Excitatory neurons_23",              #### SDK2+
                              `39`= "Excitatory neurons_11",
                              `40`= "Astrocytes_5",
                              `41`= "Inhibitory neurons_8",       #########
                              `42`= "Excitatory neurons_12",
                              `43`= "Excitatory neurons_13",
                              `44`= "Neuroendocrine cells_1",       # NeuroendocrInhibitory neuronse cells - previously identified as vSMC     ### high Inhibitory neurons Excitatory neurons markers / SDK2+ CD36+
                              `45`= "Endothelial cells_1",
                              `46`= "Inhibitory neurons_7",
                              `47`= "Oligodendrocytes precursor cells_3",
                              `48`= "Excitatory neurons_14",
                              `49`= "Excitatory neurons_15",          ######
                              `50`= "Excitatory neurons_16",
                              `51`= "Pericytes_1",
                              `52`= "Astrocytes_6",
                              `53`= "Excitatory neurons_17",
                              `54`= "Oligodendrocytes_15",
                              `55`= "Excitatory neurons_18",
                              `56`= "Excitatory neurons_19",
                              `57`= "Oligodendrocytes precursor cells_4",
                              `58`= "Oligodendrocytes precursor cells_5",
                              `59`= "Immune cells (CD45+)_3",
                              `60`= "Astrocytes_7",
                              `61`= "Oligodendrocytes_16",
                              `62`= "Oligodendrocytes precursor cells_6",
                              `63`= "Oligodendrocytes_17",
                              `64`= "Oligodendrocytes precursor cells_7",
                              `65`= "Astrocytes_9",
                              `66`= "Oligodendrocytes precursor cells_8",
                              `67`= "Excitatory neurons_20",
                              `68`= "Immune cells (CD45+)_4",
                              `69`= "Astrocytes_8",
                              `70`= "Immune cells (CD45+)_5",
                              `71`= "Astrocytes_10",
                              `72`= "Oligodendrocytes precursor cells_9",
                              `73`= "Oligodendrocytes_18",
                              `74`= "Excitatory neurons_21",
                              `75`= "Immune cells (CD45+)_6",
                              `76`= "Oligodendrocytes precursor cells_10",
                              `77`= "Excitatory neurons_22"
)

prot.combined$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"

################## IL15 expression by cell types
Idents(prot.combined) <- "cohort"
gene1 <- "IL15"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "celltypes", #cols = colors,
              idents = "Hong Kong University"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.00001, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1


FeaturePlot(prot.combined, features = gene1, raster = T, pt.size = 1.6,
            cols = c("gray80","#9E1021"),
            split.by = "cohort",
            # ncol = 2,
            # min.cutoff = 0.5,
            max.cutoff = 2.0,
            reduction = "umap") + theme(legend.position = "right")


#####################################

celltypes <- c(
  "Astrocytes",
  #"Microglia", 
  "Immune cells (CD45+)",
  "Oligodendrocytes",
  "Oligodendrocytes precursor cells",
  "Endothelial cells",
  "Pericytes",
  "Neuroendocrine cells",
  "Excitatory neurons",
  "Inhibitory neurons"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = celltypes)

col1 <- c(
  paletteer_c("ggthemes::Red-Green-Gold Diverging", 9)[9:1] 
)

DimPlot(prot.combined, reduction = "umap", raster = T, pt.size = 0.9, 
        cols = col1,
        label = T,
        repel = T,
        group.by = "celltypes",
        # split.by = "disease"
)

######## fig 5 B



myeloids <- subset(myeloids, idents = c("0","1","2","6","7",
                                             "9","10"))


DotPlot(myeloids, features = c("P2RY12","CSF1R","CD74","C3",#"CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5",
                               "IL1B","APOE","IRF8", #"CST7","CLEC7A","AIF1",
  "MS4A4A","CD163","SIGLEC1","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2","CD4", #"GDA","MARCO",
  "CD8A","CD8B","CD3E","IL2RB","KLRB1","IL7R","PRF1","GZMH","GZMA","GZMB","NKG7"#,"TRAC","TRBC1"#,"TRDC","TRGC1","GZMK","CD3D",
),
col.max = 20,
dot.scale = 10, 
cluster.idents = F, group.by = "celltypes",
#scale = F,
#split.by = "disease"
) + RotatedAxis()# + coord_flip()

myeloids <- RenameIdents(myeloids, 
                              `0`= "Microglia M0_1",
                              `1`= "Monocytes . Macrophages_1",
                              `2`= "Microglia MgND_1",
                              `6`= "Microglia M0_2",
                              `7`= "Microglia M0_3",
                              `9`= "CD8+ T cells_1",
                              `10`= "Microglia M0_4"
)
myeloids$celltypes <- sub("(.*)_.*", "\\1", Idents(myeloids))
Idents(myeloids) <- "celltypes"

celltypes <- c(
  "Microglia M0",
  "Microglia MgND",
  "Monocytes . Macrophages",
  "CD8+ T cells"
  )

myeloids$celltypes <- factor(myeloids$celltypes, levels = celltypes)

col1 <- c(
  paletteer_d("colorBlindness::SteppedSequential5Steps")[c(2,11,13,14)]
)

celltypes <- c(
  "CD8+ T cells",
  "Microglia M0",
  "Microglia MgND",
  "Monocytes / Macrophages"
)
myeloids$celltypes <- factor(myeloids$celltypes, levels = celltypes)

DimPlot(myeloids, reduction = "umap", raster = F,# pt.size = 0.9, 
        cols = col1,
        label = T,
        repel = T,
        group.by = "celltypes", 
        split.by = "cohort"
)

# saveRDS(myeloids, "/media/patrick/Portela/weiner_lab/Saef/Brain_subset_myeloids.rds")
# saveRDS(prot.combined, "/media/patrick/Portela/weiner_lab/Saef/Brain_complete_dataset.rds")
##### fig 5C

genes1 <- c("JAK1","IL15","CD36","STAT3","IRF3","TGFB1","TGFBR1",
            "IFNGR2","TLR2","MAPK14",
            "PIAS1","PPARG")
DotPlot(myeloids, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltype.disease",
        idents = c("Microglia M0_C",
                   "Microglia M0_AD",
                   "Microglia MgND_C",
                   "Microglia MgND_AD",
                   "Monocytes . Macrophages_C",
                   "Monocytes . Macrophages_AD"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

##### fig 5D

myeloids$cohort.disease <- paste(myeloids$cohort, myeloids$disease,sep = "_")
order1 <- c("UCI_C","UCI_AD","MIT_C","MIT_AD","Hong Kong University_C","Hong Kong University_AD")
myeloids$cohort.disease <- factor(myeloids$cohort.disease, levels = order1)

my_comparisons <- list(c("MIT_AD","MIT_C"),c("UCI_AD","UCI_C"),
                       c("Hong Kong University_AD","Hong Kong University_C"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1)]

p1 <- VlnPlot(myeloids, features = "IL18",
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors,
              idents = c("Monocytes / Macrophages", "Microglia MgND", "Microglia M0")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IL18"]])+0))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1

genes1 <- c("IL15")

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "cohort.disease",
                idents = "Monocytes / Macrophages"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+2))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.age2.disease.pdf"),
         plot = p1, width = 12, height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## violins

myeloids$disease <- factor(myeloids$disease, levels = c("C","AD"))
my_comparisons <- list(c("AD","C"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1)]

p1 <- VlnPlot(myeloids, features = "TGFB1",
              pt.size = 0.05, raster = F, group.by = "disease", cols = colors,
              idents = c(
                "Monocytes . Macrophages"#, 
                #"Microglia MgND"#, 
                #"Microglia M0"
                )
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["TGFB1"]])+0.5))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1


##### fig 5E

myeloids$cluster.disease <- paste(myeloids$seurat_clusters, myeloids$disease, sep = "_")
Idents(myeloids) <- "cluster.disease"

myeloids$celltype.disease <- paste(myeloids$celltypes, myeloids$disease, sep = "_")
Idents(myeloids) <- "celltype.disease"

celltypes2 <- c(
  "Microglia M0_C",
  "Microglia M0_AD",
  "Microglia MgND_C",
  "Microglia MgND_AD",
  "Monocytes . Macrophages_C",
  "Monocytes . Macrophages_AD",
  "CD8+ T cells_C",
  "CD8+ T cells_AD"
)

myeloids$celltype.disease <- factor(myeloids$celltype.disease, levels = celltypes2)

DotPlot(myeloids, features = c(#"CD8B","CD3E","IL2RB","GSDMB","KLRB1","LTB","PDCD1","IFNG","TNF",
                               "FASLG","BATF", "IFNG-AS1", #"IFNG",#"IL7R"
  "CD8A","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP" #,"IRF2","PRF1"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "celltype.disease",
idents = c("CD8+ T cells_C","CD8+ T cells_AD"),
scale = T,
#split.by = "disease"
) + RotatedAxis() + coord_flip()

DotPlot(myeloids, features = c(#"CD8B","CD3E","IL2RB","GSDMB","KLRB1","LTB","PDCD1","IFNG","TNF",
  "FASLG","BATF","IFNG-AS1", #"IL7R"
  "CD8A","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP" #,"IRF2","PRF1"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "cohort.disease",
idents = "CD8+ T cells",
scale = T,
#split.by = "cohort.disease"
) + RotatedAxis() + coord_flip()


##### Fig 5F
res <- read.xlsx("wilcox_myeloids_ADxHC_Microglia MgND_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_Microglia M0_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_CD8+ T cells_DEGs.xlsx")
res <- read.xlsx("wilcox_myeloids_ADxHC_Monocytes . Macrophages_DEGs.xlsx")
res <- read.xlsx("DEGs_pathways_ronaldo/wilcox_prot.combined_ADxHC_Immune cells (CD45+)_DEGs.xlsx")

colnames(res)[1] <- "genes"
hav <- c(#"IFNG",
  "FASLG","BATF","CD8A","GZMK","GZMH","GZMA","GZMB","NKG7","TOX","TYROBP","CX3CR1"
  # "CD8B","IL7R","PRF1","IL2RB","GSDMB","KLRB1","LTB","TNF", "PDCD1","IRF2","CD3E",
  )

hav <- c("JAK1","IL15","STAT3","TGFB1","TGFBR1","CD36","IRF3",
         "IFNGR2","TLR2","MAPK14",
         "PIAS1","PPARG"
         )

hav <- c("IFNGR2","JAK1","STAT2","STAT3","IRF1",
            "PIAS1","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP3K1","MAP3K4","MAP3K5", 
            "MAP3K14","MAP4K3",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86","IRF3","IFNGR1","MAP3K7","MAP2K5","MAPK8","MAP2K6","MAP4K4","MAP3K3","IRF2","JAK2","MAP3K20","MAPK1","MAP4K5","STAT1","RUNX3",
) # clas

EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.05,
                xlim = c(-6, 6),
                ylim = c(0, 40),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Brain myeloids - AD x C",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                maxoverlapsConnectors = Inf,#NULL
                #min.segment.length = 0.0000001,
                directionConnectors = "both",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)


##### Fig 5G

file1 <- read.delim("myeloids_up_AD_brain_pathways_sel2.csv", header = T)
colunas1 <- file1$Term[9:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3,5.5,8)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    #axis.text.x = element_text(face="bold", color = "black",size=12), 
    axis.text.y = element_text(face="bold",color = "black",size=12),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Brain myeloid cells (Microglia and monocytes) AD upregulated genes") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 



#### Figs suppl 

DotPlot(prot.combined, features = c( "GFAP","AQP4","P2RY12","CSF1R","CD74","C3",#"PTPRC",#"LCN2", "GJA1", "SLC1A2","FGFR3","NKAIN4", "SLC1A3",  #Astrocytes
                                     "MBP","MOBP","PLP1",#"MOG","CLDN11","MYRF","GALC","ERMN","MAG", "OLIG2",   #Oligodendrocytes
                                     "VCAN","PCDH15",#"CSPG4","PDGFRA", "SOX10","NEU4", "GPR37L1","C1QL1","CDO1","EPN2",   #Oligodendrocyte_precursor_cells
                                     "FLT1","CLDN5", #"VTN","ITM2A", "VWF", "FAM167B","BMX","CLEC1B",    #Endothelial_cells
                                     "PDE5A","PTH1R","P2RY14",#"AMBP","HIGD1B","COX4I2", "AOC3","ABCC9","KCNJ8","CD248",  #Pericytes
                                     "HTR2C","SLC17A7","NRGN","CAMK2A", #"SATB2", #"COL5A1","SDK2","NEFM", "SLC17A6",   #Excitatory_neurons
                                    "GAD1","GAD2"#,#"TAC1","PENK","SST","NPY","MYBPC1","PVALB","GABBR2",  "SLC32A1",  #Inhibitory_neurons
                                     #"AIF1","P2RY12","CSF1R","CD74","C3","CST3","HEXB", "C1QA", "CX3CR1","TMEM119","SLC2A5"#,   #Microglia
                                     #"CD14","FCGR3A","FCGR1A","CD68","TFRC","CCR5","ITGAM","CCR2","HP","SELL","GDA","EMILIN2", #Macrophages
                                    # "CD8A", "CD3E",
                                    
),
col.max = 20, 
dot.scale = 10, 
cluster.idents = F, group.by = "celltypes",
#scale = F,
#split.by = "cohort"
) + RotatedAxis()

FeaturePlot(myeloids, features = "CLU", raster = F, 
            ncol = 2, reduction = "umap",
            split.by = "disease",
            #max.cutoff = 3
)

FeaturePlot(prot.combined, features = "AC138207.8", raster = F, 
            ncol = 2, reduction = "umap",
            split.by = "disease",
            #max.cutoff = 3
)

FeaturePlot(prot.combined, features = c("HTR1E"), raster = F,# pt.size = 1.4, 
            # ncol = 2,
            reduction = "umap",
            # split.by = "disease",
            # max.cutoff = 2
)


order1 <- c("C","AD")
prot.combined$disease <- factor(prot.combined$disease, levels = order1)

my_comparisons <- list(c("AD","C"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1)]

p1 <- VlnPlot(prot.combined, features = "HTR1E",
              pt.size = 0.05, raster = F, group.by = "disease", cols = colors,
              idents = c("Neuroendocrine cells")
              # idents = c("Neuroendocrine cells")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][["HTR1E"]])+0.5))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.8
p1

DotPlot(prot.combined, features = c("HTR2A","HTR1E","HTR7","HTR5A-AS1","HTR2A-AS1","HTR2B"), 
        idents = "Neuroendocrine cells", cols = "RdBu",
        col.max = 20, 
        dot.scale = 10, 
        group.by = "disease") + RotatedAxis() + coord_flip()



## monocytes module score
# Idents(prot.combined) <- "orig.ident"
mono.markers <- c("MS4A4A","CD163","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2"#,"CD4","GDA","MARCO","SIGLEC1"
)

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

FeaturePlot(prot.combined, features = "Mono.markers1", raster = T, pt.size = 2.,
            cols = c("gray80","#9E1021"),
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 0.5,
            max.cutoff = 4.0,
            reduction = "umap") + theme(legend.position = "right")

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

FeaturePlot(myeloids, features = "Mono.markers1", raster = T, pt.size = 2.5,
            cols = c("gray80","#9E1021"),
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 1.5,
            max.cutoff = 3.,
            reduction = "umap") + theme(legend.position = "right")


##################################### Fig 5 - CSF

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/public_data/AD/CSF/scRNAseq_CSF_AD1/Seurat/")

prot.combined2 <- readRDS("obj_unintegrated.Rds")

prot.combined2 <- RenameIdents(prot.combined2, 
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

prot.combined2$celltypes <- sub("(.*)_.*", "\\1", Idents(prot.combined2))
Idents(prot.combined2) <- "celltypes"
prot.combined2$celltypes.disease <- paste(prot.combined2$celltypes, prot.combined2$disease, sep = "_")

DimPlot(prot.combined2, reduction = "umap.unintegrated", label = T, raster=F, repel = T
        #group.by = "condition"
        #split.by = "disease"#, ncol = 4
)

cell.types <- c(
  "Classical monocytes",
  "Intermediate monocytes",
  "mo-DC",
  "pDC",
  "Naive CD8+ T cells",
  "Effector CD8+ T cells",
  "NKT-like CD8+",
  "NK",
  "NKT-like CD4+",
  "Effector CD4+ T cells",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "B cells"
)

prot.combined2$celltypes <- factor(prot.combined2$celltypes, levels = cell.types)

DotPlot(prot.combined2, features = c(#"PF4", # platelets
                                    "CD14","FCGR3A", #monocytes
                                    "CD1C","FCER1A","IL3RA", # DCs  "ITGAX",
                                    "CD8A","CD8B", #"CD244", # T cells
                                    "PRF1","GZMM","GZMH","GZMA", # mNKs  "GZMB",
                                    "IL2RB","KLRC1","NCAM1", # NKs
                                    "KLRB1",  # NKT-like  "LILRB1", "ZBTB16",
                                    #"CD3E","CD3D", # T cells
                                    "TNF","CD4","CCR7",#"S100A4","SELL", # T cells
                                    #"IFNG", # Th1
                                    # "FOXP3", "IL2RA", # Tregs
                                    # "TRDV2","TRGC2","TRGV9", #gamma delta
                                    #"PTPRC","CCL5",
                                    #"IL7R","TBX21","EOMES",# iNKs
                                    "CD19","CD79A","CD79B","MS4A1" #,
                                    # "IGHM","IGKC","CD27","CD1D","CD22","CD86","IGLC2","IGLC3","IGHD","AIM2", "BANK1","RALGPS2","TNFRSF13B", # B cells
                                    # "IL4R","TCL1A",# Naive B cells "CXCR4", "BTG1",  "YBX3", 
                                    # "SSPN",  # Memory B cells   "COCH", "TNFRSF13C", "TEX9","LINC01781",
                                    # "LINC01857", # mature B cells
                                    # "IGHA2","TNFRSF17","TXNDC5" # plasma cells "DERL3","MZB1","POU2AF1","CPNE5","NT5DC2"
),group.by = "celltypes",
col.max = 20, 
dot.scale = 10, 
cluster.idents = F, 
) + RotatedAxis() + coord_flip()



order1 <- c("HC","MCI","AD")
prot.combined2$disease <- factor(prot.combined2$disease, levels = order1)

my_comparisons <- list(c("AD","HC"),c("AD","MCI"),c("MCI","HC"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,3,1)]

p1 <- VlnPlot(prot.combined2, features = "TOX",
              pt.size = 0.05, raster = F, group.by = "disease", cols = colors,
              idents = c("Effector CD8+ T cells")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][["TOX"]])+1.0))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.8
p1



cell.types <- c(
  # "Effector CD4+ T cells",
  # "Naive CD4+ T cells",
  "Effector CD8+ T cells",
  # "NKT-like CD4+",
  "Naive CD8+ T cells",
  "Intermediate monocytes",
  # "mo-DC",
  "Classical monocytes",
  "NK",
  # "Memory T CD4+",
  # "pDC",
  "NKT-like CD8+"#,
  # "B cells"
)

Idents(prot.combined2) <- "celltypes"
prot.combined2$celltypes.disease <- paste(prot.combined2$celltypes, prot.combined2$disease, sep = "_")
Idents(prot.combined2) <- "celltypes.disease"
### running on SCT
DefaultAssay(prot.combined2) <- "SCT"
prot.combined2 <- PrepSCTFindMarkers(prot.combined2)
for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined2, ident.1 = paste0(cell.types[i], "_AD"), 
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

### running on RNA
DefaultAssay(prot.combined2) <- "RNA"
prot.combined2 <- NormalizeData(prot.combined2, normalization.method = "LogNormalize")
# prot.combined2 <- FindVariableFeatures(prot.combined2)
# prot.combined2 <- ScaleData(prot.combined2)

for (i in 1:length(cell.types)) {
  zk.response0 <- FindMarkers(prot.combined2, ident.1 = paste0(cell.types[i], "_AD"), 
                              ident.2 = paste0(cell.types[i], "_HC"),
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
  write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_RNA_ADxHC_", cell.types[i], "_DEGs.xlsx"))
  rm(zk.response0)
}

zk.response0 <- FindMarkers(prot.combined2, ident.1 = c("NKT-like CD8+_AD", "NK_AD"),
                            #c("Intermediate monocytes_AD", "Classical monocytes_AD"), 
                            ident.2 = c("NKT-like CD8+_HC", "NK_HC"),
                            #c("Intermediate monocytes_HC", "Classical monocytes_HC"),
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
# write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_RNA_ADxHC_inter_classic.monos_DEGs.xlsx")
write.xlsx(as.data.frame(zk.response0), rowNames = T,file="wilcox_RNA_ADxHC_nkt.cd8_nks_DEGs.xlsx")
rm(zk.response0)

res <- read.xlsx("wilcox_RNA_ADxHC_Classical monocytes_DEGs.xlsx")
res <- read.xlsx("wilcox_RNA_ADxHC_Effector CD8+ T cells_DEGs.xlsx")

colnames(res)[1] <- "genes"
hav <- c(#"IFNG",
  "FASLG","BATF","CD8A","GZMK","GZMB","TOX","CX3CR1",
  "CD8B","PRF1","IL2RB","IL2RG","LTB", "PDCD1","IRF2"#,"CD3E","GSDMB","TNF","IL7R","KLRB1","NKG7","GZMH","GZMA","TYROBP",
)

hav <- c("JAK1","JAK2","IL15","IL15RA","STAT1","STAT3","TGFB1","TGFBR1","IRF1",#"IRF2",
         "IFNGR1","IFNGR2","TLR2","TLR4","MAPK14", #"CD36","PIAS1","STAT2","IRF3",
         "PPARG"#,"CREB1"#,"CREBBP"
)


EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.03,
                xlim = c(-6, 7.5),
                ylim = c(0, 12.5),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "CSF Effector CD8+ T cells - AD x C",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 200,
                maxoverlapsConnectors = Inf,#NULL
                #min.segment.length = 0.0000001,
                directionConnectors = "both",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)




genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT3","IRF1",#"TGFB1","TGFBR1",
            "TLR2","TLR4","MAPK14","IL15","IL15RA" #,
            #"PIAS1","PPARG","IRF3","CD36","STAT2","CREB1","CREBBP",
            )
DotPlot(prot.combined2, features = rev(genes1),
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltypes.disease",
        idents = c("Classical monocytes_HC",
                   #"Classical monocytes_MCI",
                   "Classical monocytes_AD"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT3","IRF1",#"TGFB1","TGFBR1",
            "TLR2","TLR4","MAPK14","IL15","IL15RA" #,
            #"PIAS1","PPARG","IRF3","CD36","STAT2","CREB1","CREBBP"
)
Idents(prot.combined2) <- "celltypes"
DotPlot(prot.combined2, features = rev(genes1),
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "disease",
        idents = c("Classical monocytes","Intermediate monocytes","mo-DC"),
        # idents = c("Intermediate monocytes_HC",
        #            "Intermediate monocytes_MCI",
        #            "Intermediate monocytes_AD"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()
Idents(prot.combined2) <- "celltypes.disease"
DotPlot(prot.combined2, features = rev(genes1),
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltypes.disease",
        idents = c("Intermediate monocytes_HC",
                   # "Intermediate monocytes_MCI",
                   "Intermediate monocytes_AD"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


FeaturePlot(prot.combined2, 
            split.by = "disease", 
            features = "IFNG", 
            raster = T,
            pt.size = 3,
            cols = c("gray80","#9E1021"),
            max.cutoff = 1.75,
            min.cutoff = 1.5
            ) + theme(legend.position = "right")

FeaturePlot(prot.combined2, 
            split.by = "disease", 
            features = "IL15", 
            raster = F,
            # pt.size = 3,
            cols = c("gray80","#9E1021"),
            # max.cutoff = 1.75,
            # min.cutoff = 1.5
) + theme(legend.position = "right")


# Idents(prot.combined) <- "orig.ident"
mono.markers <- c("MS4A4A","CD163","F13A1","MAFB","TNFAIP2","IL15","ASAH1","GAS7","PLA2G7","PLXND1","EMILIN2"#,"CD4","GDA","MARCO","SIGLEC1"
)

prot.combined2 <- AddModuleScore(
  prot.combined2,
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

FeaturePlot(prot.combined2, features = "Mono.markers1", raster = T, pt.size = 2.,
            cols = c("gray80","#9E1021"),
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 0.5,
            max.cutoff = 4.0,
            reduction = "umap.unintegrated") + theme(legend.position = "right")




###############################    Fig FTD ALS

###### Fig 3B
#####  GO bubble
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/myeloids_FTD_pathways/")

file1 <- read.delim("reactome_GOBP_FTDxPN_PFC_myeloids_upregulated_DEGs_sel.csv", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3,7,11)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    #axis.text.x = element_text(face="bold", color = "black",size=12), 
    axis.text.y = element_text(face="bold",color = "black",size=12),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Brain myeloid cells FTLD upregulated genes") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()



##################################### Fig 6 - mouse

setwd("/media/patrick/JANELSO/Bioinfo/weiner_lab/GENESIO/izzy/scRNAseq/mice/Seurat/")

prot.combined2 <- readRDS("AD-mice.izzy.no.SCT-unintegrated_annotated.rds")


## fig 6A

celltypes <- c(
  "Neutrophils",
  "Monocytes",
  "cDC1",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD8+",
  "Effector CD8+ T cells",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Cytotoxic CD4+ T cells",
  "Treg",
  "gamma-delta T cells",
  "B cells",
  "Memory B cells"
)

col1 <- c(paletteer_c("ggthemes::Green", 4),
          paletteer_c("ggthemes::Orange-Gold", 6),
          paletteer_c("ggthemes::Blue", 5),
          paletteer_c("ggthemes::Purple", 2))

prot.combined2$celltypes <- factor(prot.combined2$celltypes, levels = celltypes)

DimPlot(prot.combined2, reduction = "umap.unintegrated", raster = F,# pt.size = 0.9, 
        cols = col1,
        label = T,
        repel = T,
        group.by = "celltypes",
)


## fig 6B

DotPlot(prot.combined2, features = c(
  "S100a9","Cxcr2","Mmp9", # "S100a8","Il1r2","Retnlg", Neutrophils 
  "Cx3cr1","Fcer1g","Lyz2", #monocytes
  "Bst2", # DC "Cd1d1", "Itgax","Nrp1",
  "Ly6c1","Ifng", #   "Tnf",
                                    "Tox", #"Il2rg", iNKs 
                                    # "ITGB2","PECAM1","IL3RA","LAMP1", #baso
                                    #"Pf4", # platelets
                                    "Klrb1c","Il2rb", #"Ncam1","Klrc1",
                                    "Prf1","Gzma",#"Gzmm","Gzmk", # mNKs "Gzme",
                                    "Cd8a","Cd8b1",#"Ptprc",#"CCL5",#"CD244",
                                    "Cd4",#"S100A4","SELL", # T cells
                                    "Cd3e","Cd3d", "Ccr7",# T cells
                                    "Foxp3", "Il2ra", # Tregs
                                    "Trdc","Trgv2",#"Trgv4","Trgv5","Trgv7","Trdv1", #gamma delta
                                    "Cd19","Cd22"#,"Ighm"#,"Ighg1"#,"IGKC","IGHM","CD27","CD1D","CD22","CD86","MS4A1","IGLC2","IGLC3" #, "IGLL5" b cells
),
col.max = 20, #idents = #c("Classical Mono_1_AD","Classical Mono_2_AD","Classical Mono_1_C","Classical Mono_2_C"),
dot.scale = 10, 
cluster.idents = F, group.by = "celltypes",
#scale = F,
#split.by = "condition"
) + RotatedAxis()


####### fig 6C

order1 <- c("wt","Tg")
prot.combined2$condition <- factor(prot.combined2$condition, levels = order1)

DotPlot(prot.combined2, features = c("Il15","Il15ra","Irf2","Pias1","Tnf","Mapk14","Mapk1","Spi1","Pparg"#,   #"Creb1","Stat1","Cd36","Tlr2","Tlr4","Jak1",
                                    #"Jak2","Ifngr1","Ifngr2"
                                    #"Cd8a","Cd8b1","Irf2","Ltb","Prf1","Gzmm","Gzmk","Nkg7","Tox","Ifng","Tnf","Batf","Il2rb", "Klrb1","Il2rg"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "condition",
idents = c("Monocytes"),
scale = T,
) + RotatedAxis() + coord_flip()



####### fig 6D

DotPlot(prot.combined2, features = c(#"Il15","Il15ra","Irf2","Pias1","Tnf","Mapk14","Mapk1","Spi1","Pparg"#,   #"Creb1","Stat1","Cd36","Tlr2","Tlr4","Jak1",
                                     #"Jak2","Ifngr1","Ifngr2"
                                     "Cd8a","Cd8b1","Irf2","Ltb","Prf1","Gzmm","Gzmk","Nkg7","Tox","Ifng","Tnf","Batf","Il2rb", "Klrb1","Il2rg"
),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "condition",
idents = c("Effector CD8+ T cells"),
scale = T,
) + RotatedAxis() + coord_flip()

####### fig 6D

FeaturePlot(prot.combined2, features = "Ifng", raster = F, #split.by = "condition"
            )
FeaturePlot(prot.combined2, features = "Il15", raster = F, split.by = "condition",
            max.cutoff = 3
)


####### fig 6E

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1)]
order1 <- c("wt","Tg")
prot.combined2$condition <- factor(prot.combined2$condition, levels = order1)

my_comparisons <- list(c("wt","Tg"))
# Idents(prot.combined) <- "seurat_clusters"
p1 <- VlnPlot(prot.combined2, features = "Ifng",
              pt.size = 0.05, raster = F, group.by = "condition", cols = colors,
              idents = "iNK"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["Ifng"]])))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1




setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/mice/Seurat/")

file1 <- read.csv("pathways_up_DEGs_monocytes_Tgxwt_mouse_sel.csv", header = T)
colunas1 <- file1$Term[9:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3,5.5,8)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    #axis.text.x = element_text(face="bold", color = "black",size=12), 
    axis.text.y = element_text(face="bold",color = "black",size=12),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Mouse peripheral monocytes Tg upregulated genes") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 























