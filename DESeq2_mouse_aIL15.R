library(BiocGenerics)
library(tximport)
library(S4Vectors)
library(DESeq2)
library(openxlsx)
library(biomaRt)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(apeglm)
library(readr)
library(robustbase)
library(corrplot)
library(PerformanceAnalytics)
library(genefilter)
library(rnaseqGene)
library(EnsDb.Mmusculus.v79)
library(BiocParallel)
register(MulticoreParam(14))
library(VennDiagram)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(tibble)
library(EnhancedVolcano)
library(paletteer)

dir <- "/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/IL15_bulk/analysis/SK-6DGU/"
setwd(dir)

tx2gene <- read.delim("tx2gene.tsv", header = T)

## sample table
sampleFiles <- grep("CD8_.*_.*",list.files(dir),value=TRUE) # common info for the a specific analysis 
sampleCell <- sub("(.*)_.*_.*","\\1",sampleFiles)
samplesSex <- "male"
sampleCondition <- sub(".*_(.*)_.*","\\1",sampleFiles)
sampleReplicate <- sub(".*_.*_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          celltype = sampleCell,
                          replicate = sampleReplicate,
                          sex = samplesSex,
                          treatment = sampleCondition)
# write.csv(as.data.frame(sampleTable), file = "meta_M.csv")
# sampleTable2 <- read.csv("meta_M.csv", header = T, row.names = 1)

files <- file.path(dir, sampleTable$sampleName, "quant.sf")
names(files) <- sampleTable$sampleName

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = F)

ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, design = ~ treatment) 
  # ~ celltype)
keep <- rowSums(counts(ddsTxi)>=10) >= round(0.4*length(sampleFiles),0)
# keep <- rowSums(counts(ddsTxi)>1) >= round(0.4*length(sampleFiles),0)
ddsTxi <- ddsTxi[keep,]
# row_sub = apply(counts(ddsTxi), 1, function(row) all(row !=0 ))
# ddsTxi[row_sub,]

ddsTxi$treatment <- relevel(ddsTxi$treatment, ref = "IgG")

dds <- DESeq(ddsTxi)#, test="LRT", reduced =  ~ condition + genotype)# add ', test="LRT", reduced = ~ condition + days'  for time course, or just ', test="LRT"' for multicomparisons
res <- results(dds)
resultsNames(dds)

## to extract norm counts data
dds2 <- estimateSizeFactors(dds)
dds3 <- counts(dds2, normalized = T)

## log norm (can use 'vst' instead to center data)
vsd <- vst(dds, blind = F)
# vsd$cohort_sample <- paste(vsd$cohort, vsd$replicate, sep = "_")
# rld <- rlog(dds, blind = F)
# rld$cohort_sample <- paste(rld$cohort, rld$replicate, sep = "_")
# rld_counts <- assay(rld)

plotPCA(vsd, intgroup="treatment", ntop = 500) + ## add extra features to the plot 
  #scale_color_manual(values=c("#3333CC", "red", "#339900"))
  geom_text(aes(label=vsd$replicate), hjust=1.7) + #xlim(-23,13) + 
  ggtitle("Neurons - aIL15 vs. IgG")

## reduce cohort effect
# mat <- assay(vsd)
# mm <- model.matrix(~condition, colData(vsd))
# mat <- limma::removeBatchEffect(mat, batch=vsd$cohort, design=mm)
# assay(vsd) <- mat
# ## PCA everything
# plotPCA(vsd, intgroup="condition", ntop = 500) + ## add extra features to the plot
#   #scale_color_manual(values=c("#3333CC", "red", "#339900"))
#   geom_text(aes(label=vsd$cohort), hjust=1.1) + #xlim(-25,15) +
#   ggtitle("PDT2+3 Microglia - Females")

### Global Hierarchical clustering heatmaps (euclidean distance)
## samples distance
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds)
colnames(sampleDistMatrix) <- colnames(dds)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#png(filename = "female_microglia_global_hclus.png", units = "px", res = 300, height = 2500, width = 2500)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)
#dev.off()

## Cook's distance
par(mar=c(16,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)



setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/IL15_bulk/analysis/SK-6DGU/results/")

resultsNames(dds)
res2 <- results(dds, name = #"Intercept")
  "treatment_aIL15_vs_IgG") 
summary(res2)
res3 <- res2[!is.na(res2$padj),]
resFilter2 <- res3[res3$padj < 0.1,]
# resFilter2 <- resFilter2[resFilter2$baseMean > 20,]
# resFilter2 <- resFilter2[abs(resFilter2$log2FoldChange) > 13,]
markers2 <- dds3[rownames(dds3) %in% rownames(resFilter2),]
# targets <- c(
#   #"S100a9","Cxcr2","Mmp9", # "S100a8","Il1r2","Retnlg", #Neutrophils
#   "Cx3cr1","Fcer1g", #monocytes "Lyz2",
#   "Cd1d1",# DC "Nrp1", "Itgax", "Bst2",
#   "Ifng", "Tnf",# "Ly6c1",
#   "Tox", "Il2rg", #iNKs
#   # "ITGB2","PECAM1","IL3RA","LAMP1", #baso
#   #"Pf4", # platelets
#   "P2ry12","Csf1r", "Cst3", "Hexb", "C1qa", "Cx3cr1", "Aif1", "Tmem119",  #microglia "C3", "Iba1","Cd74",
#   "Slc17a6",  "Slc17a7",  "Nrgn","Nefm",    #excitatory_neurons "Sdk2", "Col5a1","Camk2a","Satb2",
#   "Slc32a1",  "Gad1", "Gad2", "Tac1", "Penk", "Sst",  "Npy",  "Mybpc1", "Pvalb", "Gabbr2",   #inhibitory_neurons
#   "Klrb1c","Il2rb", "Ncam1","Klrc1",
#   "Prf1","Gzma","Gzmm","Gzmk","Gzme","Gzmb", # mNKs
#   "Cd8a","Cd8b1","Ptprc",#"CCL5",#"CD244",
#   #"Cd4",#"S100A4","SELL", # T cells
#   "Cd3e","Cd3d","Ccr7", # T cells
#   #"Foxp3", "Il2ra", # Tregs
#   "Trdc","Trgv2","Trgv4","Trgv5","Trgv7","Trdv1"#, #gamma delta
#   #"Cd19","Cd22","Ighm","Ighg1"#,"IGKC","IGHM","CD27","CD1D","CD22","CD86","MS4A1","IGLC2","IGLC3" #, "IGLL5" b cells
# )
# markers2 <- dds3[rownames(dds3) %in% targets,]

sampleTable2 <- data.frame(
  #replicate = sampleTable$replicate,
  treatment = sampleTable$treatment,
  # celltype = sampleTable$celltype,
  # cohort = sampleTable$cohort,
  row.names = sampleTable$sampleName)

## Targets heatmap
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdGy")) )(paletteLength)
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
heat <- pheatmap(markers2, cluster_rows = T, show_rownames = T, 
                 show_colnames = F, annotation_col = sampleTable2,
                 cluster_cols = F, annotation_names_col = T, annotation_legend = T,
                 color = colors,
                 border_color = NA,
                 main = "Neurons",
                 #labels_row = as.expression(newnames),
                 angle_col = 45,
                 scale = "row"
) 
#dev.off()


### extract stats
res2Ordered <- res2[order(res2$padj),]
#resFilter2Ordered <- resFilter2[order(resFilter2$log2FoldChange),]
# write.xlsx(as.data.frame(markers2), rowNames = T, file = "PDT7_mic_females_WT_Naivex3KL_A_DEGs_norm.counts.xlsx")
write.xlsx(as.data.frame(dds3), rowNames = T, file = "norm.counts_mouse_PS19_Neurons_ail15xigg.xlsx")
write.xlsx(as.data.frame(res2Ordered), rowNames = T, file="mouse_PS19_Neurons_ail15xigg_stats.xlsx")
#write.xlsx(as.data.frame(resFilter2Ordered), rowNames = T, file="PDT2_males_CxB_DEGs_stats.xlsx")

write.xlsx(as.data.frame(dds3), rowNames = T, file = "basemean.filter.norm.counts_mouse_PS19_Neurons_ail15xigg.xlsx")
write.xlsx(as.data.frame(res2Ordered), rowNames = T, file="basemean.filter.mouse_PS19_Neurons_ail15xigg_stats.xlsx")





###### heatmap highlighting only target genes

colors <- colorRampPalette( paletteer_c("grDevices::Tropic", 30)  )(paletteLength)
# colors <- colorRampPalette( rev(paletteer_c("grDevices::Blue-Yellow 3", 30) ) )(paletteLength)
# colors <- colorRampPalette( rev(paletteer_c("grDevices::Green-Orange", 30) ) )(paletteLength)
# colors <- colorRampPalette( rev(paletteer_c("grDevices::ArmyRose", 30) ) )(paletteLength)
# colors <- colorRampPalette( rev(paletteer_c("grDevices::Earth", 30)) )(paletteLength)
# colors <- colorRampPalette( rev(paletteer_c("grDevices::TealRose", 30)) )(paletteLength)
# colors <- colorRampPalette( paletteer_c("viridis::cividis", 30) )(paletteLength)
# colors <- colorRampPalette( paletteer_c("grDevices::Cividis", 30) )(paletteLength)
heat <- pheatmap(markers2, cluster_rows = T, show_rownames = T, 
                 show_colnames = F, annotation_col = sampleTable2,
                 cluster_cols = F, annotation_names_col = T, annotation_legend = T,
                 color = colors,
                 border_color = NA,
                 main = "CD8+ T cells - aIL15 vs. IgG",
                 #labels_row = as.expression(newnames),
                 # angle_col = 45,
                 scale = "row"
) 
#dev.off()

add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}

# hav2 <- c("H2-Aa","H2-Ab1","H2-DMa","H2-DMb2","H2-Eb1","H2-Ob","Il4ra","Cd83","Ccr7","Cd24a","Dusp10","Trav19") # B
hav2 <- c("Stat3","Prf1","Gzmb","Gzma","Tnfaip3","Klrb1b","Klrb1c","Ifngr1","Il18rap","Irf8","Klri2") # CD8



hav2 <- t(read.delim("microglia_targets.txt", header = F)) 
hav2 <- t(read.delim("neuron_targets2.txt", header = F))

add.flag(heat, kept.labels = hav2, repel.degree = 0.5)
#### Volcano 1
EnhancedVolcano(res2Ordered,
                lab = row.names(res2Ordered), #NA,
                selectLab = hav2,#tmarkers[1,],
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.000983653608845477,
                FCcutoff = 0.5,
                xlim = c(-6, 7),
                ylim = c(0, 12),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 6.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Neurons - aIL15 vs. IgG",
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




#####################           Volcanos 2
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/IL15_bulk/analysis/SK-6DGU/results/new_names/")
# hav <- read.xlsx("mouse_PS19_cd8_ail15xigg_stats.xlsx")
# hav <- read.xlsx("mouse_PS19_microglia_ail15xigg_stats.xlsx")
hav <- read.xlsx("mouse_PS19_Neurons_ail15xigg_stats.xlsx")

hav <- hav[!is.na(hav$padj),]
colnames(hav)[1] <- "genes"
# selected_geness <- c("Stat3","Prf1","Gzmb","Gzma","Tnfaip3","Klrb1b","Klrb1c","Ifngr1","Il18rap","Irf8","Klri2") # CD8
# selected_geness <- t(read.delim("microglia_targets.txt", header = F))
selected_geness <- t(read.delim("neuron_targets.txt", header = F))

hav$neglog10p <- -log10(hav$padj)

# Assign color categories based on your thhavholds
hav$color <- "Not significant"
hav$color[hav$padj < 0.1 & hav$log2FoldChange > 0] <- "Upregulated"
hav$color[hav$padj < 0.1 & hav$log2FoldChange < 0] <- "Downregulated"

# Get coordinates for selected geness
label_data <- hav[hav$genes %in% selected_geness, ]

# Split label data into up- and down-regulated
up <- label_data[label_data$log2FoldChange > 0, ]
down <- label_data[label_data$log2FoldChange < 0, ]

# Define stacking positions
x_label_right <- max(hav$log2FoldChange) - 1.4
x_label_left <- min(hav$log2FoldChange) + 1.5
y_range <- range(hav$neglog10p)
if (nrow(up) > 0) {
  y_label_up <- seq(y_range[2], y_range[1], length.out = nrow(up))
  up$label_x <- x_label_right
  up$label_y <- y_label_up
}
if (nrow(down) > 0) {
  y_label_down <- seq(y_range[2], y_range[1], length.out = nrow(down))
  down$label_x <- x_label_left
  down$label_y <- y_label_down
}
stack_labels <- rbind(up, down)

# Volcano plot with color by threshold and stacked labels
ggplot(hav, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = color), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "salmon", "Not significant" = "grey60", "Downregulated" = "lightblue")) +
  # Connectors
  geom_segment(
    data = stack_labels,
    aes(x = log2FoldChange, y = neglog10p, xend = label_x, yend = label_y),
    color = "gray40", linewidth = 0.7, inherit.aes = FALSE
  ) +
  # Stacked labels
  geom_label(
    data = stack_labels,
    aes(x = label_x, y = label_y, label = genes),
    hjust = ifelse(stack_labels$log2FoldChange > 0, 0, 1),
    fontface = "bold", fill = "white", color = "black", label.size = 0.3, inherit.aes = FALSE
  ) +
  xlab(expression("Log"[2]*"(Fold-Change)")) +
  ylab(expression("-Log"[10]*"("*italic("p")*"-adjusted)")) +
  ggtitle("Neurons - aIL15 vs. Untreated") +
  coord_cartesian(xlim = c(x_label_left - 1, x_label_right + 1)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())



########## pathways

file1 <- read.delim("mouse_PS19_cd8_ail15xigg_GOBP_down_select.txt", header = T)
file1 <- read.delim("mouse_PS19_mic_ail15xigg_pathways_down_select.txt", header = T)
file1 <- read.delim("mouse_PS19_neurons_ail15xigg_pathways_up_select.txt", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3.3,4.5,5.7)
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
  ggtitle("PS19 Neurons - aIL15 vs. Untreated (Upregulated DEGs)") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()








#############################   Nichenet   #########################
library(nichenetr)
library(tidyverse)

 
## sample table Neuron
sampleFiles <- grep("Neuron_.*_.*",list.files(dir),value=TRUE) # common info for the a specific analysis 
sampleCell <- sub("(.*)_.*_.*","\\1",sampleFiles)
samplesSex <- "male"
sampleCondition <- sub(".*_(.*)_.*","\\1",sampleFiles)
sampleReplicate <- sub(".*_.*_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          celltype = sampleCell,
                          replicate = sampleReplicate,
                          sex = samplesSex,
                          treatment = sampleCondition)
# write.csv(as.data.frame(sampleTable), file = "meta_M.csv")
# sampleTable2 <- read.csv("meta_M.csv", header = T, row.names = 1)

files <- file.path(dir, sampleTable$sampleName, "quant.sf")
names(files) <- sampleTable$sampleName

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = F)

ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, design = ~ treatment)
# keep <- rowSums(counts(ddsTxi)>=10) >= round(0.4*length(sampleFiles),0)
keep <- rowSums(counts(ddsTxi)>1) >= round(0.4*length(sampleFiles),0)
ddsTxi <- ddsTxi[keep,]
# row_sub = apply(counts(ddsTxi), 1, function(row) all(row !=0 ))
# ddsTxi[row_sub,]

ddsTxi$treatment <- relevel(ddsTxi$treatment, ref = "IgG")

dds_Neuron <- DESeq(ddsTxi)#, test="LRT", reduced =  ~ condition + genotype)# add ', test="LRT", reduced = ~ condition + days'  for time course, or just ', test="LRT"' for multicomparisons
res_Neuron <- results(dds_Neuron)
# resultsNames(dds)

## to extract norm counts data
dds2_Neuron <- estimateSizeFactors(dds_Neuron)
dds3_Neuron <- counts(dds2_Neuron, normalized = T)

resultsNames(dds_Neuron)
res2 <- results(dds_Neuron, name = "treatment_aIL15_vs_IgG") 
summary(res2)
res3 <- res2[!is.na(res2$padj),]
resFilter2 <- res3[res3$padj < 0.1,]
# resFilter2 <- resFilter2[resFilter2$baseMean > 20,]
neurons_up <- row.names(resFilter2[resFilter2$log2FoldChange > 0,])
neurons_down <- row.names(resFilter2[resFilter2$log2FoldChange < 0,])

## sample table CD8
sampleFiles <- grep("^CD8_.*_.*",list.files(dir),value=TRUE) # common info for the a specific analysis 
sampleCell <- sub("(.*)_.*_.*","\\1",sampleFiles)
samplesSex <- "male"
sampleCondition <- sub(".*_(.*)_.*","\\1",sampleFiles)
sampleReplicate <- sub(".*_.*_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          celltype = sampleCell,
                          replicate = sampleReplicate,
                          sex = samplesSex,
                          treatment = sampleCondition)
# write.csv(as.data.frame(sampleTable), file = "meta_M.csv")
# sampleTable2 <- read.csv("meta_M.csv", header = T, row.names = 1)

files <- file.path(dir, sampleTable$sampleName, "quant.sf")
names(files) <- sampleTable$sampleName

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = F)

ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, design = ~ treatment)
keep <- rowSums(counts(ddsTxi)>=10) >= round(0.4*length(sampleFiles),0)
# keep <- rowSums(counts(ddsTxi)>1) >= round(0.4*length(sampleFiles),0)
ddsTxi <- ddsTxi[keep,]
# row_sub = apply(counts(ddsTxi), 1, function(row) all(row !=0 ))
# ddsTxi[row_sub,]

ddsTxi$treatment <- relevel(ddsTxi$treatment, ref = "IgG")

dds_CD8 <- DESeq(ddsTxi)#, test="LRT", reduced =  ~ condition + genotype)# add ', test="LRT", reduced = ~ condition + days'  for time course, or just ', test="LRT"' for multicomparisons
res_CD8 <- results(dds_CD8)
# resultsNames(dds)

## to extract norm counts data
dds2_CD8 <- estimateSizeFactors(dds_CD8)
dds3_CD8 <- counts(dds2_CD8, normalized = T)

resultsNames(dds_CD8)
res2 <- results(dds_CD8, name = "treatment_aIL15_vs_IgG") 
summary(res2)
res3 <- res2[!is.na(res2$padj),]
resFilter2 <- res3[res3$padj < 0.1,]
# resFilter2 <- resFilter2[resFilter2$baseMean > 20,]
cd8_up <- row.names(resFilter2[resFilter2$log2FoldChange > 0,])
cd8_down <- row.names(resFilter2[resFilter2$log2FoldChange < 0,])


## sample table microglia
sampleFiles <- grep("^FCRLS_.*_.*",list.files(dir),value=TRUE) # common info for the a specific analysis 
sampleCell <- sub("(.*)_.*_.*","\\1",sampleFiles)
samplesSex <- "male"
sampleCondition <- sub(".*_(.*)_.*","\\1",sampleFiles)
sampleReplicate <- sub(".*_.*_(.*)","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          celltype = sampleCell,
                          replicate = sampleReplicate,
                          sex = samplesSex,
                          treatment = sampleCondition)
# write.csv(as.data.frame(sampleTable), file = "meta_M.csv")
# sampleTable2 <- read.csv("meta_M.csv", header = T, row.names = 1)

files <- file.path(dir, sampleTable$sampleName, "quant.sf")
names(files) <- sampleTable$sampleName

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = F)

ddsTxi <- DESeqDataSetFromTximport(txi.salmon, colData = sampleTable, design = ~ treatment)
keep <- rowSums(counts(ddsTxi)>=10) >= round(0.4*length(sampleFiles),0)
# keep <- rowSums(counts(ddsTxi)>1) >= round(0.4*length(sampleFiles),0)
ddsTxi <- ddsTxi[keep,]
# row_sub = apply(counts(ddsTxi), 1, function(row) all(row !=0 ))
# ddsTxi[row_sub,]

ddsTxi$treatment <- relevel(ddsTxi$treatment, ref = "IgG")

dds_mic <- DESeq(ddsTxi)#, test="LRT", reduced =  ~ condition + genotype)# add ', test="LRT", reduced = ~ condition + days'  for time course, or just ', test="LRT"' for multicomparisons
res_mic <- results(dds_mic)
# resultsNames(dds)

## to extract norm counts data
dds2_mic <- estimateSizeFactors(dds_mic)
dds3_mic <- counts(dds2_mic, normalized = T)

resultsNames(dds_mic)
res2 <- results(dds_mic, name = "treatment_aIL15_vs_IgG") 
summary(res2)
res3 <- res2[!is.na(res2$padj),]
resFilter2 <- res3[res3$padj < 0.1,]
# resFilter2 <- resFilter2[resFilter2$baseMean > 20,]
mic_up <- row.names(resFilter2[resFilter2$log2FoldChange > 0,])
mic_down <- row.names(resFilter2[resFilter2$log2FoldChange < 0,])



# 3. Prepare NicheNet input sets ----------------------------------
### start analysis
organism = "mouse"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
}
# ligand_target_matrix <- nichenet_ligand_target_matrix
ligand_receptor_network <- lr_network

# Predict active ligands
ligands <- unique(ligand_receptor_network$from)






# For interaction: one cell type as "sender", other as "receiver"
# CD8+ T cells are sender, Neurons are receiver
expressed_genes_sender <- rownames(dds3_CD8)[rowMeans(dds3_CD8) > 1]
expressed_genes_receiver <- rownames(dds3_Neuron)[rowMeans(dds3_Neuron) > 1]

geneset_oi <- neurons_down # genes downregulated in neurons with treatment

potential_ligands <- intersect(ligands, cd8_down)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
write.xlsx(as.data.frame(ligand_activities), rowNames = T, file="CD8_neurons_nichenet_ligands.xlsx")

# View top predicted ligands
head(ligand_activities)

best_upstream_ligands = ligand_activities %>% top_n(121, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()


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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="CD8_2_neuron_nichenet_vis_ligand_target.xlsx")



file1 <- ligand_activities
colnames(file1)[1] <- "ligands"
# colunas1 <- file1$ligands[c(4,5,2,3,6,1)]
colunas1 <- file1$ligands[c(1,6,3,2,5,4)]
mid <- 0
# file1$logperson <- -log10(file1$pearson)
ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(ligands, levels = colunas1), 
                                      size = auroc, fill = aupr_corrected), alpha = 1, shape = 21) +
  scale_size(range = c(3, 8), name = expression("AUROC"), #breaks = c(20,60,100)
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
  ggtitle("CD8+ T cells -> Neurons (downregulated CCC)") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,1)* max(abs(file1$aupr_corrected)), 
                       name = expression("corrected AUPR"#"Log"[2]*"(Fold-Enrichment)"
                       )) 





markers2 <- dds3_CD8[rownames(dds3_CD8) %in% colunas1,]
markers2 <- markers2[colunas1[6:1],]
sampleTable2 <- data.frame(
  #replicate = sampleTable$replicate,
  treatment = sampleTable$treatment,
  # celltype = sampleTable$celltype,
  # cohort = sampleTable$cohort,
  row.names = sampleTable$sampleName)

## Targets heatmap
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
heat <- pheatmap(markers2, cluster_rows = F, show_rownames = T, 
                 show_colnames = F, annotation_col = sampleTable2,
                 cluster_cols = F, annotation_names_col = T, annotation_legend = T,
                 color = colors,
                 border_color = NA,
                 main = "CD8+ T cells ligands",
                 #labels_row = as.expression(newnames),
                 angle_col = 45,
                 scale = "row"
) 
#dev.off()




# Microglia are sender, Neurons are receiver
expressed_genes_sender <- rownames(dds3_mic)[rowMeans(dds3_mic) > 1]
expressed_genes_receiver <- rownames(dds3_Neuron)[rowMeans(dds3_Neuron) > 1]

geneset_oi <- neurons_down # genes downregulated in neurons with treatment

potential_ligands <- intersect(ligands, mic_down)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
write.xlsx(as.data.frame(ligand_activities), rowNames = T, file="mic_2_neurons_nichenet_ligands2.xlsx")

# View top predicted ligands
head(ligand_activities)

best_upstream_ligands = ligand_activities %>% top_n(121, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()


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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="mic_2_neuron_nichenet_vis_ligand_target2.xlsx")



file1 <- ligand_activities
# file1 <- file1[c(1:7,9:11,13:18),]
colnames(file1)[1] <- "ligands"
file1 <- file1[order(file1$aupr_corrected),]
colunas1 <- file1$ligands[c(1:16)]
# colunas1 <- file1$ligands[c(11,15,18,2,4,17,14,10,13,7,9,6,3,1,5,16)]
mid <- 0
# file1$logperson <- -log10(file1$pearson)
ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(ligands, levels = colunas1), 
                                      size = auroc, fill = aupr_corrected), alpha = 1, shape = 21) +
  scale_size(range = c(3, 8), name = expression("AUROC"), #breaks = c(20,60,100)
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
  ggtitle("Microglia -> Neurons (downregulated CCC)") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,1)* max(abs(file1$aupr_corrected)), 
                       name = expression("corrected AUPR"#"Log"[2]*"(Fold-Enrichment)"
                       )) 









# Microglia are sender, CD8 are receiver
expressed_genes_sender <- rownames(dds3_mic)[rowMeans(dds3_mic) > 1]
expressed_genes_receiver <- rownames(dds3_CD8)[rowMeans(dds3_CD8) > 1]

geneset_oi <- cd8_down # genes downregulated in CD8 with treatment

potential_ligands <- intersect(ligands, mic_down)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
write.xlsx(as.data.frame(ligand_activities), rowNames = T, file="mic_2_cd8_nichenet_ligands2.xlsx")

# View top predicted ligands
head(ligand_activities)

best_upstream_ligands = ligand_activities %>% top_n(121, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()


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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="mic_2_cd8_nichenet_vis_ligand_target2.xlsx")



file1 <- ligand_activities
file1 <- file1[c(1:7,9:11,13:18),]
colnames(file1)[1] <- "ligands"
file1 <- file1[order(file1$aupr_corrected),]
colunas1 <- file1$ligands[c(1:18)]
# colunas1 <- file1$ligands[c(11,15,18,2,4,17,14,10,13,7,9,6,3,1,5,16)]
mid <- 0
# file1$logperson <- -log10(file1$pearson)
ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(ligands, levels = colunas1), 
                                      size = auroc, fill = aupr_corrected), alpha = 1, shape = 21) +
  scale_size(range = c(3, 8), name = expression("AUROC"), #breaks = c(20,60,100)
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
  ggtitle("Microglia -> CD8+ T cells (downregulated CCC)") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,1)* max(abs(file1$aupr_corrected)), 
                       name = expression("corrected AUPR"#"Log"[2]*"(Fold-Enrichment)"
                       )) 









# CD8+ T cells are sender, microglia are receiver
expressed_genes_sender <- rownames(dds3_CD8)[rowMeans(dds3_CD8) > 1]
expressed_genes_receiver <- rownames(dds3_mic)[rowMeans(dds3_mic) > 1]

geneset_oi <- mic_down # genes downregulated in neurons with treatment

potential_ligands <- intersect(ligands, cd8_down)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
write.xlsx(as.data.frame(ligand_activities), rowNames = T, file="CD8_mic_nichenet_ligands.xlsx")

# View top predicted ligands
head(ligand_activities)

best_upstream_ligands = ligand_activities %>% top_n(121, aupr_corrected) %>% 
  arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()


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
write.xlsx(as.data.frame(vis_ligand_target2), rowNames = T, file="CD8_2_mic_nichenet_vis_ligand_target.xlsx")



file1 <- ligand_activities
colnames(file1)[1] <- "ligands"
file1 <- file1[order(file1$aupr_corrected),]
colunas1 <- file1$ligands[c(1:18)]
# colunas1 <- file1$ligands[c(4,5,2,3,6,1)]
# colunas1 <- file1$ligands[c(1,6,3,2,5,4)]
mid <- 0
# file1$logperson <- -log10(file1$pearson)
ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(ligands, levels = colunas1), 
                                      size = auroc, fill = aupr_corrected), alpha = 1, shape = 21) +
  scale_size(range = c(3, 8), name = expression("AUROC"), #breaks = c(20,60,100)
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
  ggtitle("CD8+ T cells -> Microglia (downregulated CCC)") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,1)* max(abs(file1$aupr_corrected)), 
                       name = expression("corrected AUPR"#"Log"[2]*"(Fold-Enrichment)"
                       )) 











# ## bulksignalR
# library(BulkSignalR)
# library(tidyverse)
# library(biomaRt)
# library(mygene)
# 
# # 1. Load bulk RNA-seq data -------------------------------------
# # Assume input files are genes x samples, columns: neurons_control1, ..., cd8_treat3, etc.
# neurons <- dds3_Neuron[row.names(dds3_Neuron) %in% row.names(dds3_CD8),]
# cd8 <- dds3_CD8[row.names(dds3_CD8) %in% row.names(neurons),]
# 
# # Combine counts for all samples (or TPM/FPKM if normalized; counts preferred)
# expr_data_mouse <- cbind(neurons, cd8)
# 
# 
# # ------------------------
# # 2. Map mouse genes to human gene symbols using mygene
# # ------------------------
# genes_mouse <- row.names(expr_data_mouse)
# 
# # Directly use queryMany() without creating MyGeneInfo object
# query_res <- queryMany(genes_mouse,
#                        scopes = "symbol",
#                        fields = c("homologene"),  # to get homolog info
#                        species = "mouse",
#                        # return.as = "DataFrame",
#                        returnall = TRUE)
# 
# # Extract human homologene symbols
# # Each result may contain homologene entries; extract human ortholog symbol (taxid 9606)
# extract_human_symbol <- function(hit) {
#   if (!is.null(hit$homologene$genes)) {
#     homologs <- hit$homologene$genes
#     # homologs is a list of vectors: [[taxid, gene_symbol]]
#     human_symbol <- NA_character_
#     for (gene in homologs) {
#       if (gene[1] == 9606) {
#         human_symbol <- gene[2]
#         break
#       }
#     }
#     return(human_symbol)
#   } else {
#     return(NA_character_)
#   }
# }
# 
# human_symbols <- sapply(query_res$hits, extract_human_symbol)
# mouse_symbols <- sapply(query_res$hits, function(x) if(!is.null(x$query)) x$query else NA_character_)
# 
# # Create mapping dataframe and remove NA human symbols
# mapping_df <- data.frame(mouse=mouse_symbols, human=human_symbols, stringsAsFactors = FALSE)
# mapping_df <- mapping_df[!is.na(mapping_df$human) & mapping_df$human != "", ]
# 
# # Filter expression data: keep only mouse genes with mapped human orthologs
# expr_data_mouse_filtered <- expr_data_mouse[rownames(expr_data_mouse) %in% mapping_df$mouse, ]
# 
# # Rename rownames using mapped human symbols
# human_names_for_mouse <- mapping_df$human[match(rownames(expr_data_mouse_filtered), mapping_df$mouse)]
# rownames(expr_data_mouse_filtered) <- human_names_for_mouse
# 
# # Optionally, remove duplicated human gene symbols by averaging rows
# expr_data_human <- expr_data_mouse_filtered %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "human_symbol") %>%
#   group_by(human_symbol) %>%
#   summarize(across(.cols = everything(), .fns = mean)) %>%
#   column_to_rownames(var = "human_symbol") %>%
#   as.matrix()
# 
# # ------------------------
# # 3. Prepare metadata
# # ------------------------
# sample_names <- colnames(expr_data_human)
# celltype <- sub("(.*)_.*_.*","\\1",sample_names)
# condition <- sub(".*_(.*)_.*","\\1",sample_names)
# metadata <- data.frame(sample=sample_names, celltype=celltype, condition=condition, stringsAsFactors = FALSE)
# 
# 
# bsrdm <- prepareDataset(counts = expr_data_human)
# 
# 
# 
# 
# 
# # 4. Create the BulkSignalR Data Model object
# bsrdm <- BSRDataModel(counts = expr_data_human,
#                                    species = "hsapiens",        # because data is human-symbol mapped
#                                    # If dealing with mouse original genes, supply conversion.dict & species="mmusculus"
#                                    # conversion.dict = conversion_dict_df,
#                                    # species = "mmusculus",
#                                    colData = metadata)
# 
# # Learn model parameters (fast mode for testing)
# bsrdm <- BulkSignalR::learnParameters(bsrdm, quick = TRUE)
# 
# # Run ligand-receptor inference using Reactome database
# bsrinf <- BulkSignalR::BSRInference(bsrdm, reference = "REACTOME")
# 
# # Optionally reset gene names to initial organism using your conversion dictionary if available:
# # bsrinf <- BulkSignalR::resetToInitialOrganism(bsrinf, conversion.dict = conversion_dict_df)
# 
# # Visualize results
# BulkSignalR::plotSignatureHeatmap(bsrinf)
# 
# # Explore inference object bsrinf for top ligand-receptor interactions
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # ------------------------
# # 4. Run BulkSignalR analysis
# # ------------------------
# signal <- BulkSignalR::Signal(expr_data_human, metadata, group.col="celltype", condition.col="condition")
# 
# # Sender: CD8 T cells, Receiver: Neurons
# lr_cd8_to_neuron <- BulkSignalR::findLR(signal,
#                                         sender="CD8",
#                                         receiver="Neuron",
#                                         condition.test="aIL15",
#                                         condition.ref="IgG",
#                                         only.diff=TRUE,
#                                         pval.cutoff=0.05,
#                                         logfc.cutoff=0.5)
# 
# # Sender: Neurons, Receiver: CD8 T cells
# lr_neuron_to_cd8 <- BulkSignalR::findLR(signal,
#                                         sender="Neuron",
#                                         receiver="CD8",
#                                         condition.test="aIL15",
#                                         condition.ref="IgG",
#                                         only.diff=TRUE,
#                                         pval.cutoff=0.05,
#                                         logfc.cutoff=0.5)
# 
# # ------------------------
# # 5. Visualize results
# # ------------------------
# BulkSignalR::plotEnrichedLR(lr_cd8_to_neuron)
# BulkSignalR::plotEnrichedLR(lr_neuron_to_cd8)
# 
# # ------------------------
# # 6. Export results and show top interactions
# # ------------------------
# write.csv(lr_cd8_to_neuron, "BulkSignalR_cd8_to_neuron_results.csv", row.names=FALSE)
# write.csv(lr_neuron_to_cd8, "BulkSignalR_neuron_to_cd8_results.csv", row.names=FALSE)
# 
# head(lr_cd8_to_neuron)
# head(lr_neuron_to_cd8)



































