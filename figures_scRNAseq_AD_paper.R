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


setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/")

prot.combined <- readRDS("./obj_harmony_patient_reg.out.mito.ncounts.v3c.Rds")

########### monocytes not different

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

prot.combined$cluster.disease <- paste(prot.combined$seurat_clusters, prot.combined$disease, sep = "_")

prot.combined$cluster.disease <- factor(prot.combined$cluster.disease, levels = c(
  "50_C","50_AD","47_C","47_AD","29_C","29_AD","7_C","7_AD",
  "68_C","68_AD","71_C","71_AD","43_C","43_AD","44_C","44_AD"
))

Idents(prot.combined) <- "cluster.disease"

IL15.sig <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
              "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
              "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
              "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
              "CREB1","CREBBP",
              "IL15","IL15RA"#,"CD86"
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

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1)]

gene1 <- "IL15"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cluster.disease", cols = colors,
              idents = #classic
                c(
                  "50_C","50_AD","47_C","47_AD","29_C","29_AD","7_C","7_AD",
                  "68_C","68_AD","71_C","71_AD","43_C","43_AD","44_C","44_AD"
                )
              # c("50","47","29","7","68","71","43","44")
                # "Classical monocytes"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1










################# CD8 status
order1 <- c("male_C","male_AD","female_C","female_AD")
prot.combined$sex.disease <- factor(prot.combined$sex.disease, levels = order1)
my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)
DotPlot(prot.combined, features = c("CXCL10","CXCR3","CXCL16","CXCR6"),
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, 
        group.by = "sex.disease",
        # group.by = "age.range.disease2",
        idents = c("Effector CD8+ T cells"),
        scale = T,
        # scale.by = "size",
        #split.by = "disease"
) + RotatedAxis() + coord_flip()

FeaturePlot(prot.combined, raster = F,# pt.size = 1.5,
            features = c("CXCL16"), 
            # features = c("CXCR6","CXCL16"),
            split.by = "disease",
            # min.cutoff = 1,
            # max.cutoff = 2.5,
            cols = c("gray90","#9E1021"),
            reduction = "umap") + theme(legend.position = "right")

gene1 <- "CXCR3"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "disease", #cols = colors,
              idents = "Effector CD8+ T cells"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1


####################       APOE status
IL15.sig <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
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
              pt.size = 0.05, raster = F, group.by = "geno.disease2", #cols = colors,
              idents = "Classical monocytes"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1

cytotox <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
  "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
  "NKG7","TYROBP",#"FASLG","KLRB1",
  "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
  "ENTPD1"
  #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
) # effector cd8

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









meta.data <- read.csv("meta.txt", header = T)
prot.combined$geno2 <- prot.combined$patient

for (i in 1:length(meta.data$GEM)) {
  prot.combined$geno2 <- recode(prot.combined$geno2, "meta.data$GEM[i] = meta.data$APOE.geno2[i]")
}
prot.combined$disease.geno2 <- paste(prot.combined$disease, prot.combined$geno2, sep = "_")
prot.combined$sex.geno2 <- paste(prot.combined$sex, prot.combined$geno2, sep = "_")
prot.combined$age.geno2 <- paste(prot.combined$age.range2, prot.combined$geno2, sep = "_")
prot.combined <- subset(prot.combined, subset = geno2 %in% c("High","Neut"))
prot.combined <- subset(prot.combined, subset = disease %in% c("AD"))

my_comparisons <- c("C_Reduced","C_Neut","C_High","AD_Reduced","AD_Neut","AD_High")
prot.combined$disease.geno2 <- factor(prot.combined$disease.geno2, levels = my_comparisons)
prot.combined$age.geno2 <- factor(prot.combined$age.geno2, levels = c("50-70_Neut","50-70_High","71-90_Neut","71-90_High"))
prot.combined$sex.geno2 <- factor(prot.combined$sex.geno2, levels = c("male_Neut","male_High","female_Neut","female_High"))

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 10,
        # cluster.idents = T,
        #scale = F,
        group.by = "sex.geno2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
  "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
  "NKG7","TYROBP",#"FASLG","KLRB1",
  "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
  "ENTPD1"
  #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
  
) # effector cd8

DotPlot(prot.combined, features = rev(genes1),
col.max = 20,
cols = "RdBu",
dot.scale = 10, 
cluster.idents = F, group.by = "age.geno2",
idents = c("Effector CD8+ T cells"),
scale = T,
# scale.by = "size",
#split.by = "disease"
) + RotatedAxis() + coord_flip()

gene1 <- "IL15"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "patient", #cols = colors,
              idents = "Classical monocytes"
                # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1


######################## corr IL15 x Tau and IL15 x ab42
######## corr plots

library(ggpubr)
mydf1 <- read.xlsx("Il15_Tau_Corelation_6-3-25_control.xlsx",rowNames = F)

ggscatter(mydf1, x = "IL15", y = "Tau", #label = rownames(mydf1),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "pearson", label.x = 5, label.y = 140) +
  annotate("text", x = 10, y = 150, label = (paste0("slope==", coef(lm(mydf1$Tau~mydf1$IL15))[2])), parse = T) + 
  ggtitle(expression("IL15 x Tau - Healthy controls")) + ylab(expression("Tau")) +
  xlab(expression("IL15"))

ggscatter(mydf1, x = "IL15", y = "Tau", #label = rownames(mydf1),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "pearson", label.x = 20, label.y = 230) +
  annotate("text", x = 27, y = 240, label = (paste0("slope==", coef(lm(mydf1$Tau~mydf1$IL15))[2])), parse = T) + 
  ggtitle(expression("IL15 x Tau - AD")) + ylab(expression("Tau")) +
  xlab(expression("IL15"))



mydf1 <- read.xlsx("ab42_IL15_ELISA_human_7-11-25_control.xlsx",rowNames = F)

ggscatter(mydf1, x = "IL15", y = "Ab42", #label = rownames(mydf1),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "pearson", label.x = 5, label.y = 160) +
  annotate("text", x = 10, y = 170, label = (paste0("slope==", coef(lm(mydf1$Ab42~mydf1$IL15))[2])), parse = T) + 
  ggtitle(expression("IL15 x Ab42 - Healthy controls")) + ylab(expression("Ab42")) +
  xlab(expression("IL15"))

mydf1 <- read.xlsx("ab42_IL15_ELISA_human_7-11-25_AD.xlsx",rowNames = F)

ggscatter(mydf1, x = "IL15", y = "Ab42", #label = rownames(mydf1),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "black", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
  stat_cor(method = "pearson", label.x = 20, label.y = 360) +
  annotate("text", x = 27, y = 370, label = (paste0("slope==", coef(lm(mydf1$Ab42~mydf1$IL15))[2])), parse = T) + 
  ggtitle(expression("IL15 x Ab42 - AD")) + ylab(expression("Ab42")) +
  xlab(expression("IL15"))


###############################        fig 1     ##############################

######## Fig 1A
# scatter plot sex and age

meta.data <- read.csv("meta.txt", header = T)
males <- meta.data[meta.data$sex == "male",]
males$y <- runif(39, 0, 1)

ggplot(males, aes(x = Age, y = y, color = Disease)) +
  geom_text(label = "\u2642", family = "Arial Unicode MS", size = 5) +
  scale_x_continuous(
    limits = c(50, 90),
    breaks = seq(50, 90, by = 10),
  ) +
  #scale_color_manual(values = c("red", "blue"))+
  scale_color_brewer(palette="Dark2") + 
  labs(title = "",
       x = "Age (years)",
       y = "Male") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # remove y-axis labels
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_blank())  # remove y-axis ticks

females <- meta.data[meta.data$sex == "female",]
females$y <- runif(31, 0, 1)

ggplot(females, aes(x = Age, y = y, color = Disease)) +
  geom_text(label = "\u2640", family = "Arial Unicode MS", size = 5) +
  scale_x_continuous(
    limits = c(50, 90),
    breaks = seq(50, 90, by = 10),
  ) +
  scale_color_brewer(palette="Dark2") + 
  labs(title = "",
       x = "Age (years)",
       y = "Female") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),  # remove y-axis labels
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks.y = element_blank())  # remove y-axis ticks

# stacked bar plot

# myCol <- paletteer_c("grDevices::Dynamic", 14) 
# myCol <- paletteer_d("ggsci::category20c_d3", 14) 
# myCol <- paletteer_d("ggsci::category20b_d3",14) 
# myCol <- paletteer_d("ggsci::default_igv",14) 
# meta.data <- read.csv("meta2.csv", header = T)
# meta.data <- meta.data %>% mutate(across(everything(), as.character))
# df_tidy <- meta.data %>% pivot_longer(-GEM, names_to = "ID")
# df_tidy$count <- rep(1, 350)
# 
# 
# p <- ggplot(data = df_tidy, aes(x=ID,y=count, fill= value, label = value) ) + 
#   geom_col(position = "stack") + scale_fill_manual(values = myCol) +
#   theme_light(base_size = 18) + 
#   #geom_text(aes(label = value),  position = position_stack(0.5)) + 
#   coord_flip()
# p

# # alluvial
# meta.data <- read.csv("meta2.csv", header = T)
# #meta.data$count <- rep(1/70, 70)
# ord <- list(meta.data$n3, NULL, NULL,NULL,NULL)
# alluvial(
#   select(meta.data, 
#          Disease, sex, Age, genotype, cohort),
#   freq=meta.data$n1,
#   col = meta.data$colors, #ifelse(test$EXP.Level == "DOWN", "azure4", "lightgray"),
#   #border = test$colors,#ifelse(test$EXP.Level == "DOWN", "azure4", "lightgray"),
#   layer = meta.data$Disease,
#   alpha = 0.7,
#   blocks="bookends",
#   gap.width = 0.06,
#   cw = 0.16,
#   cex.axis = 1.2,
#   cex = 1,
#   xw = 0.15,
#   ordering = ord,
#   axis_labels = c("Disease", "sex", "Age", "Genotype","cohort")
# ) 

# histogram
meta.data <- read.csv("meta.txt", header = T)
males <- meta.data[meta.data$sex == "female",]
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1)]

ggplot(males, aes(x = Age)) +
  geom_histogram(aes(y = ..density..), binwidth = 4, fill = "skyblue", alpha = 0.5, position = "identity", color = "black") +
  geom_density(aes(color = Disease), linewidth = 1) +
  scale_color_manual(values = c("C" = colors[2], "AD" = colors[1])) +
  scale_x_continuous(
    limits = c(50, 90),
    breaks = seq(50, 90, by = 10),
  ) +
  labs(title = "Females",
       x = "Age (years)",
       y = "Density") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1)
  )



# Create the histogram
ggplot(males, aes(x = Age, fill = Disease)) +
  geom_histogram(aes(y = ..density..), binwidth = 10, position = "dodge", 
                 alpha = 0.7) +
  scale_fill_manual(values = c("AD" = colors[2], "C" = colors[1])) +
  scale_x_continuous(
    limits = c(50, 90),
    breaks = seq(50, 90, by = 10),
  ) +
  labs(title = "Males",
       x = "Age (years)",
       y = "Density",
       fill = "Disease") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1)
  )


# Create the histogram
ggplot(males, aes(x = Age, fill = Disease)) +
  geom_histogram(aes(y = ..density..), binwidth = 5, position = "identity", 
                 alpha = 0.5) +
  scale_fill_manual(values = c("AD" = colors[2], "C" = colors[1])) +
  scale_x_continuous(
    limits = c(50, 90),
    breaks = seq(50, 90, by = 10),
  ) +
  labs(title = "Males",
       x = "Age (years)",
       y = "Density",
       fill = "Sex") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1)
  )

# Create the overlapped histograms
ggplot(males, aes(x = Age, fill = Disease)) +
  geom_histogram(aes(y = ..density..), binwidth = 4, position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("AD" = colors[2], "C" = colors[1])) +
  labs(title = "Males",
       x = "Age (years)",
       y = "Density",
       fill = "Disease") +
  theme_minimal() +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black", linewidth = 1)
  )



####### Fig 1B

myCol <- brewer.pal(3, "Set2")

DimPlot(prot.combined, reduction = "umap.unintegrated", raster = F, 
        label = F, repel = T,
        cols = myCol,
        group.by = "cohort",
) #, combine = F)


###### Fig 1C

clusters <- as.character(t(read.delim("cluster.order.txt", header = F)))

prot.combined$seurat_clusters <- factor(prot.combined$seurat_clusters, levels = clusters)

col2 <- c(
  paletteer_c("ggthemes::Green-Gold", 24),
  paletteer_c("grDevices::Heat", 27), 
  paletteer_c("grDevices::TealGrn", 19), 
  paletteer_c("grDevices::PurpOr", 13) )

# highlight monocytes
mono  <- c("50","47","29","7","68","71","43","44","19","58","18","48","36") # "69" has 1 cell in HC 80 and 78 absent in HC
mono  <- c("50","47","29","7","68","71","43","69","44","80","19","58","78","18","48","36") 
col1 <- c(paletteer_c("ggthemes::Green-Gold", 24))
col2 <- c(rep("grey50",4),col1[5:11],"grey50",col1[13],"grey50",col1[15:19],rep("grey50",64)) # without 69, 80 and 78
col2 <- c(rep("grey50",4),col1[5:20],rep("grey50",63))

# highlight nks nkts 
NK.NKT <- c("1","2","67","57","21","41","62","20","3","30","52","72","12","28","17","25","64")
col1 <- c(paletteer_c("grDevices::Heat", 27))
col2 <- c(rep("grey50",24),col1[1:8],rep("grey50",3),col1[12:20],rep("grey50",39))

# highlight effects cd8 
effect.cd8 <- c("9","14")
col1 <- c(paletteer_c("grDevices::Heat", 27))
col2 <- c(rep("grey50",44),col1[c(14,22)],rep("grey50",37))

DimPlot(prot.combined, reduction = "umap", raster = F,# pt.size = 1.3, 
        # cols = col2,
        label = T,
        repel = T,
        group.by = "seurat_clusters",
)

###### Fig 1D

## Broad cell markers

celltypes3 <- c("Monocytes","Platelets","NKT-like","NK","gamma-delta T cells",
                "CD4+ T cells","CD8+ T cells","B cells","DC")
Idents(prot.combined) <- "celltypes3"
prot.combined$celltypes3 <- factor(prot.combined$celltypes3, levels = celltypes3)

DotPlot(prot.combined, features = c(
  "CD14","FCGR3A","PF4","GZMH","NCAM1","IL2RB","KLRB1","CD4","CCR7","CD8A","CD8B",
  "IGHM","CD22","CD1C","FCER1A"
),group.by = "celltypes3",
col.max = 20, 
dot.scale = 10, 
cluster.idents = F, 
) + RotatedAxis()

## Specific cell markers
celltypes <- c(
  "Platelets",
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD4+",
  "NKT-like CD8+",
  "Effector CD8+ T cells",
  "Exhausted CD8+ T cells",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Th1",
  "Cytotoxic CD4+ T cells",
  "Treg",
  "gamma-delta T cells",
  "Naive B cells",
  "Memory B cells",
  "Cytotoxic B cells",
  "mature B cells",
  "Plasma B cells"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = celltypes)

DotPlot(prot.combined, features = c("PF4", # platelets
                                    "CD14","FCGR3A", #monocytes
                                    "ITGAX","CD1C","FCER1A","IL3RA", # DCs
                                    "KLRC1","NCAM1","IL2RB", # NKs
                                    "PRF1","GZMM","GZMH","GZMA","GZMB", # mNKs
                                    "LILRB1", "KLRB1", "ZBTB16", # NKT-like
                                    #"CD3E","CD3D", # T cells
                                    "CD8A","CD8B", #"CD244", # T cells
                                    "CD4","CCR7",#"S100A4","SELL", # T cells
                                    "TNF",#"IFNG", # Th1
                                    "FOXP3", "IL2RA", # Tregs
                                    "TRDV2","TRGC2","TRGV9", #gamma delta
                                   #"PTPRC","CCL5",
                                   #"IL7R","TBX21","EOMES",# iNKs
                                    "IGHM",#"CD19","IGKC","CD27","CD1D","CD22","CD86","MS4A1","IGLC2","IGLC3","IGHD","CD79A","CD79B","AIM2", "BANK1","RALGPS2","TNFRSF13B", # B cells
                                    "IL4R","TCL1A",# Naive B cells "CXCR4", "BTG1",  "YBX3", 
                                    "SSPN", "TEX9","LINC01781", # Memory B cells   "COCH", "TNFRSF13C", 
                                    "LINC01857", # mature B cells
                                    "IGHA2","TNFRSF17","TXNDC5" # plasma cells "DERL3","MZB1","POU2AF1","CPNE5","NT5DC2"
),group.by = "celltypes",
col.max = 20, 
dot.scale = 10, 
cluster.idents = F, 
) + RotatedAxis()




###### Fig 1E

celltypes <- c(
  "Platelets",
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD4+",
  "NKT-like CD8+",
  "Effector CD8+ T cells",
  "Exhausted CD8+ T cells",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Th1",
  "Cytotoxic CD4+ T cells",
  "Treg",
  "gamma-delta T cells",
  "Naive B cells",
  "Memory B cells",
  "Cytotoxic B cells",
  "mature B cells",
  "Plasma B cells"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = celltypes)

col1 <- c(paletteer_c("ggthemes::Green", 6),
          paletteer_c("ggthemes::Orange-Gold", 8),
          paletteer_c("ggthemes::Blue", 6),
          paletteer_c("ggthemes::Purple", 5))

DimPlot(prot.combined, reduction = "umap", raster = T, pt.size = 1.5, 
        cols = col1,
        label = F,
        repel = T,
        group.by = "celltypes",
        split.by = "cohort"
) + theme(legend.position = "none")

###### Fig 1F

# DEGs by cell type
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
new_row <- data.frame(clusters="Exhausted CD8+ T cells", DEGs=0,upregulated=0,downregulated=0)
new_df <- rbind(data[1:13, ], new_row, data[(14:24), ])
coluna1 <- new_df$clusters
coluna1 <- rev(coluna1) 
new_df$clusters <- factor(new_df$clusters, levels=coluna1)


axis_margin <- 5.5

p1 <- ggplot(new_df, aes(y=factor(clusters,levels=coluna1), x=downregulated)) +
  geom_col(fill = "slategray1") +
  geom_text(aes(label = downregulated), hjust = 1.5) + 
  scale_x_reverse(
    limits = c(2200, 0),
    breaks = seq(1500, 0, by = -500)
  ) + 
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
    axis.text.y = element_blank(),
  )

p2 <- ggplot(new_df, aes(y=factor(clusters,levels=coluna1), x=upregulated)) +
  geom_col(fill = "salmon") +
  geom_text(aes(label = upregulated), hjust = -0.5) + 
  scale_x_continuous(
    #limits = c(0, 1500),
    #breaks = seq(0, 1500, by = 500),
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
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
  )

ggarrange(p1, p2)

# ## by sex / age
# 
# data <- data.frame(matrix(ncol = 4,nrow = 0))
# 
# for (i in 1:length(celltypes)){
#   list1 <- read.xlsx(paste0("DEGs/age2.sex.cohort.regress/celltypes_age2/MAST_71-90_ADxHC_",celltypes[i],"_DEGs.xlsx"))
#   data[i,1] <- celltypes[i]
#   data[i,2] <- length(row.names(list1))
#   data[i,3] <- length(list1$avg_log2FC[list1$avg_log2FC > 0])
#   data[i,4] <- length(list1$avg_log2FC[list1$avg_log2FC < 0])
#   rm(list1)
# }
# 
# colnames(data) <- c("clusters","DEGs","upregulated","downregulated")
# 
# coluna1 <- data$clusters
# coluna1 <- rev(coluna1) 
# data$clusters <- factor(data$clusters, levels=coluna1)
# 
# axis_margin <- 5.5
# 
# p1 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x=upregulated)) +
#   geom_col(fill = "darkorange") +
#   scale_x_reverse() + 
#   scale_y_discrete(position = "right") +
#   theme(
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
#   ) + theme(
#     ## Set background color to white
#     panel.background = element_rect(fill = "white"),
#     ## Set the color and the width of the grid lines for the horizontal axis
#     panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
#     ## Remove tick marks by setting their length to 0
#     axis.ticks.length = unit(0, "mm"),
#     axis.text.y = element_blank(),
#   )
# 
# p2 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x=downregulated)) +
#   geom_col(fill = "skyblue") + scale_x_continuous(
#     #limits = c(0, 1500),
#     #breaks = seq(0, 1500, by = 500),
#   ) +
#   theme(
#     axis.title.y = element_blank(),
#     plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
#     axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
#   ) + theme(
#     ## Set background color to white
#     panel.background = element_rect(fill = "white"),
#     ## Set the color and the width of the grid lines for the horizontal axis
#     panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
#     ## Remove tick marks by setting their length to 0
#     axis.ticks.length = unit(0, "mm"),
#     ## Remove labels from the vertical axis
#     axis.text.y = element_blank(),
#   )
# 
# ggarrange(p1, p2)
# 
# 
# 
# 
# #### sex in same plot
# 
# data <- data.frame(matrix(ncol = 4,nrow = 0))
# 
# for (i in 1:length(celltypes)){
#   list1 <- read.xlsx(paste0("DEGs/age2.sex.cohort.regress/celltypes_sex/MAST_females_ADxHC_",celltypes[i],"_DEGs.xlsx"))
#   data[i,1] <- celltypes[i]
#   data[i,2] <- length(row.names(list1))
#   data[i,3] <- length(list1$avg_log2FC[list1$avg_log2FC > 0])
#   data[i,4] <- length(list1$avg_log2FC[list1$avg_log2FC < 0])
#   rm(list1)
# }
# 
# colnames(data) <- c("clusters","DEGs","upregulated","downregulated")
# 
# coluna1 <- data$clusters
# coluna1 <- rev(coluna1) 
# data$clusters <- factor(data$clusters, levels=coluna1)
# 
# data <- data[,c(1,3,4)]
# data <- data %>% pivot_longer(cols = c(upregulated, downregulated), names_to = "Variable", values_to = "Value")
# 
# data2 <- data.frame(matrix(ncol = 4,nrow = 0))
# 
# for (i in 1:length(celltypes)){
#   list1 <- read.xlsx(paste0("DEGs/age2.sex.cohort.regress/celltypes_sex/MAST_males_ADxHC_",celltypes[i],"_DEGs.xlsx"))
#   data2[i,1] <- celltypes[i]
#   data2[i,2] <- length(row.names(list1))
#   data2[i,3] <- length(list1$avg_log2FC[list1$avg_log2FC > 0])
#   data2[i,4] <- length(list1$avg_log2FC[list1$avg_log2FC < 0])
#   rm(list1)
# }
# 
# colnames(data2) <- c("clusters","DEGs","upregulated","downregulated")
# 
# coluna1 <- data2$clusters
# coluna1 <- rev(coluna1) 
# data2$clusters <- factor(data2$clusters, levels=coluna1)
# 
# data2 <- data2[,c(1,3,4)]
# data2 <- data2 %>% pivot_longer(cols = c(upregulated, downregulated), names_to = "Variable", values_to = "Value")
# 
# axis_margin <- 1.5
# myCol <- rep(c("darkorange","skyblue"),24)
# 
# p1 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x = Value, fill = Variable)) +
#   geom_col(fill = myCol) +
#   scale_x_reverse(limits = c(1800, 0),
#                   breaks = seq(1250, 0, by = -250)) + 
#   scale_y_discrete(position = "right") +
#   theme(
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     plot.margin = margin(axis_margin, 0, axis_margin, axis_margin)
#   ) + theme(
#     ## Set background color to white
#     panel.background = element_rect(fill = "white"),
#     ## Set the color and the width of the grid lines for the horizontal axis
#     panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
#     ## Remove tick marks by setting their length to 0
#     axis.ticks.length = unit(0, "mm"),
#     axis.text.y = element_blank()
#   ) + xlab("DEGs - females")
# 
# p2 <- ggplot(data2, aes(y=factor(clusters,levels=coluna1), x = Value, fill = Variable)) +
#   geom_col(fill = myCol) + scale_x_continuous(
#     limits = c(0, 1400),
#     breaks = seq(0, 1250, by = 250),
#   ) +
#   theme(
#     axis.title.y = element_blank(),
#     plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
#     axis.text.y.left = element_text(margin = margin(0, axis_margin, 0, axis_margin))
#   ) + theme(
#     ## Set background color to white
#     panel.background = element_rect(fill = "white"),
#     ## Set the color and the width of the grid lines for the horizontal axis
#     panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
#     ## Remove tick marks by setting their length to 0
#     axis.ticks.length = unit(0, "mm"),
#     ## Remove labels from the vertical axis
#     axis.text.y = element_blank()
#   ) + xlab("DEGs - males")
# ggarrange(p1, p2)



###### Fig 1G and 1H

celltypes <- c(
  "Platelets",
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "NKs and NKTs",
  # "iNK",
  # "mNK",
  # "NKT-like CD4+",
  # "NKT-like CD8+",
  "CD8+ T cells",
  # "Effector CD8+ T cells",
  # "Exhausted CD8+ T cells",
  # "Naive CD8+ T cells",
  # "Memory T CD8+",
  "CD4+ T cells",
  # "Naive CD4+ T cells",
  # "Memory T CD4+",
  # "Th1",
  # "Cytotoxic CD4+ T cells",
  # "Treg",
  "gamma-delta T cells",
  "B cells",
  # "Naive B cells",
  # "Memory B cells",
  # "Cytotoxic B cells",
  # "mature B cells",
  "Plasma B cells"
)

col1 <- c(paletteer_c("ggthemes::Green", 6),
          paletteer_c("ggthemes::Orange-Gold", 2), # paletteer_c("ggthemes::Orange-Gold", 8),
          paletteer_c("ggthemes::Blue", 2), #paletteer_c("ggthemes::Blue", 6),
          paletteer_c("ggthemes::Purple", 2) #paletteer_c("ggthemes::Purple", 5)
          )

data <- read.delim("Summary_ligand_receptor_counts_AD_samples_cellphoneDB_2.tsv", header = T)
data$n_interactions <- data$n_interactions / sum(data$n_interactions)
data$n_interactions <- data$n_interactions**3.5 
colnames(data) <- c("source_ligand","target_receptor","interaction_count")
data$source_ligand <- factor(data$source_ligand, levels = celltypes)

colnames(data) <- c("source_ligand","target_receptor","interaction_count")

adj_matrix <- xtabs(interaction_count ~ source_ligand + target_receptor, data = data)

# cols <- colorRamp2(range(adj_matrix), c("#F8F9FF", "#0E3F5C"))
cols <- colorRamp2(range(adj_matrix), c("white", "black")) #"darkred"
# cols <- colorRamp2(c(0.0,9.6e-08), c("white","darkred"))

# Reset circos parameters after plotting
circos.clear()
# Reduce plot size
par(mar = c(1, 1, 1, 1))  # Reduce margins
circos.par(
  "canvas.xlim" = c(-1.6, 1.6),  # Reduce x-axis limits
  "canvas.ylim" = c(-1.6, 1.6)   # Reduce y-axis limits
)


# Create chord diagram
chordDiagram(
  adj_matrix,
  col = cols,
  grid.col = col1,
  transparency = 0.2,
  link.lwd = 1,
  link.lty = 1,
  annotationTrack = "grid"
  #link.border = "white"
)

# Add node labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(-0.1, 0.5)
    )
  },
  bg.border = NA
)

# # Add title
title("AD cell-cell interactions")

# Add color legend
lgd = Legend(col_fun = cols, 
             title = "Interaction rate",
             # at = seq(min(adj_matrix), max(adj_matrix), length.out = 5))
             at = round(seq(min(adj_matrix), max(adj_matrix), length.out = 5),9))

draw(lgd, x = unit(1, "npc") - unit(2, "cm"), y = unit(1, "npc") - unit(2, "cm"))

# Reset circos parameters after plotting
circos.clear()





## interaction rate fold-change targeted

data <- read.delim("Summary_ligand_receptor_counts_AD_HC_samples_cellphoneDB_target_3.tsv", header = T)
data <- data[data$source_ligand %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"),]
data <- data[data$target_receptor %in% c("iNK","mNK","NKT-like CD8+","Effector CD8+ T cells"),]
# data$n_interactions <- data$n_interactions / sum(data$n_interactions)
data <- data[c("source_ligand","target_receptor","FC.HC")]

colnames(data) <- c("source_ligand","target_receptor","interaction_count")

adj_matrix <- xtabs(interaction_count ~ source_ligand + target_receptor, data = data)

# cols <- colorRamp2(range(adj_matrix), c("#F8F9FF", "#0E3F5C"))
# cols <- colorRamp2(range(adj_matrix), c("lightblue","lightpink")) #"darkred"
cols <- colorRamp2(c(0.8,0.95,1,1.05,1.2,1.45,1.7), c("lightblue","#D9E6F7","white","#FFC2C5","lightpink","indianred1","darkred"))

# Reset circos parameters after plotting
circos.clear()
# Reduce plot size
par(mar = c(1, 1, 1, 1))  # Reduce margins
circos.par(
  "canvas.xlim" = c(-1.6, 1.6),  # Reduce x-axis limits
  "canvas.ylim" = c(-1.6, 1.6)   # Reduce y-axis limits
)

col1 <- c(paletteer_c("ggthemes::Green", 3),
          paletteer_c("ggthemes::Orange-Gold", 4) # paletteer_c("ggthemes::Orange-Gold", 8),
          # paletteer_c("ggthemes::Blue", 2), #paletteer_c("ggthemes::Blue", 6),
          # paletteer_c("ggthemes::Purple", 2) #paletteer_c("ggthemes::Purple", 5)
)


# Create chord diagram
chordDiagram(
  adj_matrix,
  col = cols,
  grid.col = col1,
  transparency = 0.,
  link.lwd = 1,
  link.lty = 1,
  annotationTrack = "grid"
  #link.border = "white"
)

# Add node labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(-0.1, 0.5)
    )
  },
  bg.border = NA
)

# # Add title
title("Cell-cell interactions")

# Add color legend
lgd = Legend(col_fun = cols, 
             title = "Interaction rate fold-change",
             # at = seq(min(adj_matrix), max(adj_matrix), length.out = 5))
             at = round(seq(min(adj_matrix), max(adj_matrix), length.out = 5),1))

draw(lgd, x = unit(0.9, "npc") - unit(2, "cm"), y = unit(1, "npc") - unit(2, "cm"))

# Reset circos parameters after plotting
circos.clear()




## interaction rate by group
data <- read.delim("Summary_ligand_receptor_counts_AD_HC_samples_cellphoneDB_target_3.tsv", header = T)
data <- data[data$source_ligand %in% c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes"),]
data <- data[data$target_receptor %in% c("iNK","mNK","NKT-like CD8+","Effector CD8+ T cells"),]
# data$n_interactions <- data$n_interactions / sum(data$n_interactions)
data <- data[c("source_ligand","target_receptor","HC.percent")]

colnames(data) <- c("source_ligand","target_receptor","interaction_count")

adj_matrix <- xtabs(interaction_count ~ source_ligand + target_receptor, data = data)

# cols <- colorRamp2(range(adj_matrix), c("white","darkred")) #"darkred"
# cols <- colorRamp2(c(0.01,0.09), c("white","darkred"))
cols <- colorRamp2(c(0,0.021,0.025,0.05,0.09), c("white","#FFE7E8","lightpink2","indianred1","darkred"))
# cols <- colorRamp2(c(0,0.021,0.03,0.09), c("white","#FFE7E8","indianred1","darkred"))

# Reset circos parameters after plotting
circos.clear()
# Reduce plot size
par(mar = c(1, 1, 1, 1))  # Reduce margins
circos.par(
  "canvas.xlim" = c(-1.6, 1.6),  # Reduce x-axis limits
  "canvas.ylim" = c(-1.6, 1.6)   # Reduce y-axis limits
)

col1 <- c(paletteer_c("ggthemes::Green", 3),
          paletteer_c("ggthemes::Orange-Gold", 4) # paletteer_c("ggthemes::Orange-Gold", 8),
          # paletteer_c("ggthemes::Blue", 2), #paletteer_c("ggthemes::Blue", 6),
          # paletteer_c("ggthemes::Purple", 2) #paletteer_c("ggthemes::Purple", 5)
)


# Create chord diagram
chordDiagram(
  adj_matrix,
  col = cols,
  grid.col = col1,
  transparency = 0.,
  link.lwd = 1,
  link.lty = 1,
  annotationTrack = "grid"
  #link.border = "white"
)

# Add node labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(-0.1, 0.5)
    )
  },
  bg.border = NA
)

# # Add title
title("HC cell-cell interactions")

# Add color legend
lgd = Legend(col_fun = cols, 
             title = "Relative interaction rate",
             # at = seq(min(adj_matrix), max(adj_matrix), length.out = 5))
             # at = round(seq(min(adj_matrix), max(adj_matrix), length.out = 5),3))
             at = seq(0.01, 0.09, length.out = 5))

draw(lgd, x = unit(0.9, "npc") - unit(2, "cm"), y = unit(1, "npc") - unit(2, "cm"))

# Reset circos parameters after plotting
circos.clear()


###############################        fig 2     ##############################

###### Fig 2A
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



# myeloids <- subset(prot.combined, idents = c("Classical monocytes","Intermediate monocytes","Nonclassical monocytes","mo-DC"))
Idents(prot.combined) <- "seurat_clusters"
myeloids <- subset(prot.combined, subset = seurat_clusters %in% effect.cd8)
Idents(myeloids) <- "seurat_clusters"

myeloids$seurat_clusters <- factor(myeloids$seurat_clusters, levels = effect.cd8)

col2 <- c(#paletteer_c("ggthemes::Green", 8),
          #paletteer_c("grDevices::TealGrn", 2),
          #paletteer_c("grDevices::TealGrn", 5)
  #paletteer_c("ggthemes::Green", 13),
  #paletteer_c("ggthemes::Green-Gold", 13)#, # mono
  paletteer_c("ggthemes::Orange-Gold", 2)#, NK.NKT
  #paletteer_c("grDevices::Heat", 2)#, 
  #paletteer_c("ggthemes::Blue", 19),
  #paletteer_c("grDevices::TealGrn", 19), 
  #paletteer_c("ggthemes::Purple", 13)
  #paletteer_c("grDevices::PurpOr", 13) 
  )

DimPlot(myeloids, reduction = "umap", 
        raster = T,
        pt.size = 1.3,
        cols = col2,
        label = T,
        repel = T,
        group.by = "seurat_clusters"
)



###### Fig 2B

# my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
# prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)

Idents(prot.combined) <- "celltypes"
genes1 <- c("IFNG","TNF","TGFB1","HMGB1") 

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        #group.by = "age.range.disease2",
        #idents = NK.NKT,
) + RotatedAxis() #+ coord_flip()

celltypes <- c(
  "Platelets",
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD4+",
  "NKT-like CD8+",
  "Effector CD8+ T cells",
  "Exhausted CD8+ T cells",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Th1",
  "Cytotoxic CD4+ T cells",
  "Treg",
  "gamma-delta T cells",
  "Naive B cells",
  "Memory B cells",
  "Cytotoxic B cells",
  "mature B cells",
  "Plasma B cells"
)

prot.combined$celltypes <- factor(prot.combined$celltypes, levels = celltypes)

p1 <- VlnPlot(prot.combined, features = "HMGB1", pt.size = 0.01, raster = T,
        group.by = "celltypes"#, ncol = 2#stack = T, flip = T
        ) + theme(legend.position = "none")
# p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
#   scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IFNG"]])))+
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

FeaturePlot(prot.combined, raster = T, pt.size = 3, 
            features = genes1, 
            split.by = "disease",
            # min.cutoff = 0,
            max.cutoff = 3,
            reduction = "umap") #+ theme(legend.position = "right")
cols1 <- paletteer_c("grDevices::Reds 3", 500)[500:1] 

FeaturePlot(prot.combined, raster = T, pt.size = 1.5,
            features = "IL15", 
            split.by = "disease",
            # min.cutoff = 1,
            max.cutoff = 2.0,
            reduction = "umap") + theme(legend.position = "right")

FeaturePlot(prot.combined, raster = F, #T, pt.size = 1.5,
            features = "IFNG", 
            split.by = "disease",
            # min.cutoff = 1,
            max.cutoff = 2.5,
            cols = c("gray80","#9E1021"),
            reduction = "umap") + theme(legend.position = "right")

FeaturePlot(prot.combined, raster = T, pt.size = 1.5,
            features = "IL15", 
            # split.by = "disease",
            # min.cutoff = 1,
            # max.cutoff = 2.5,
            cols = c("gray80","#9E1021"),
            reduction = "umap") + theme(legend.position = "right")

Idents(prot.combined) <- "seurat_clusters"

prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")
order1 <- c("50-70_male_C","50-70_male_AD",
            "50-70_female_C","50-70_female_AD",
            "71-90_male_C","71-90_male_AD",
            "71-90_female_C","71-90_female_AD")
prot.combined$sex.age2.disease <- factor(prot.combined$sex.age2.disease, levels = order1)

my_comparisons <- list(c("50-70_male_AD","50-70_male_C"),c("71-90_male_AD","71-90_male_C"),
                       c("50-70_female_AD","50-70_female_C"),c("71-90_female_AD","71-90_female_C"))

my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1)]
my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
p1 <- VlnPlot(prot.combined, features = "IFNG", pt.size = 0.01, raster = F,
        # ncol = 2, stack = T, flip = T
        group.by = "sex.age2.disease", cols = colors,
        idents = c("52")#,"41,"52","62"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IFNG"]])+1.5)) +
 stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1)]
gene1 <- "IFNG"
my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
p1 <- VlnPlot(prot.combined, features = gene1, pt.size = 0.01, raster = F,
              # ncol = 2, stack = T, flip = T
              group.by = "sex.age2.disease", cols = colors,
              idents = c("52")#,"41,"52","62"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[gene1]])+1.5)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

Idents(prot.combined) <- "seurat_clusters"
col2 <- c(rep("grey50",41),"#FFFF00",rep("grey50",10),"#FFFF00",rep("grey50",9),"#FFFF00",rep("grey50",20))
col2 <- c("#FFFF00",rep("grey50",52),rep("grey50",30))
order1 <- as.character(c(52,0:51,53:82))
# order1 <- as.character(c(0:82))
prot.combined$seurat_clusters <- factor(prot.combined$seurat_clusters, levels = order1)

DimPlot(prot.combined, reduction = "umap", raster = T, pt.size = 2,
        cols = col2,
        label = F,
        repel = T,
        group.by = "seurat_clusters",
) + theme(legend.position = "none")

Idents(prot.combined) <- "seurat_clusters"
cells <- CellsByIdentities(prot.combined, idents = c(52))
DimPlot(prot.combined, reduction = "umap", raster = F, #pt.size = 1,
        cells.highlight = cells[setdiff(names(cells), "NA")],
        cols.highlight = "#FFFF00",
        #cols = col2,
        label = F,
        repel = T,
        sizes.highlight = 0.01,
        # alpha = ,
        # group.by = "seurat_clusters",
        split.by = "disease"
) + theme(legend.position = "none")






# my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
# prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)
# 
# colors <- brewer.pal(n=3,name = "Dark2")
# colors <- colors[c(2,1,2,1)]
# my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
# p1 <- VlnPlot(prot.combined, features = "TSPO", pt.size = 0.01, raster = F,
#               # ncol = 2, stack = T, flip = T
#               group.by = "age.range.disease2", cols = colors,
#               idents = "Classical monocytes"
# ) + theme(legend.position = "none")
# p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
#   scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["TSPO"]])+1.0)) +
#   stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
# p1$layers[[2]]$aes_params$alpha <- 0.1
# p1
# 





###### Fig 2A and 2B

data <- read.delim("AD_group_CellChat_cluster_level_communication_IFNG_IFNGR1_IFNGR2_only.tsv", header = T)
# data.1 <- unique(data$source)
# write.csv(data.1, file = "clusters_1.csv")
# clus.cells0 <- read.delim("cell_clusters.txt", header = F)
# clus.cells <- clus.cells0$V1
# data$source <- factor(data$source, levels = data.1)

col2 <- c(paletteer_c("ggthemes::Green-Gold", 24)[1:16],
          paletteer_c("grDevices::Heat", 27)[1:8], 
          paletteer_c("grDevices::TealGrn", 23)[1:2], 
          paletteer_c("grDevices::PurpOr", 13)[1:12],
          paletteer_c("ggthemes::Green-Gold", 24)[17:24],
          paletteer_c("grDevices::Heat", 27)[9:17],
          paletteer_c("grDevices::TealGrn", 23)[3:5],
          paletteer_c("grDevices::Heat", 27)[18:24],
          paletteer_c("grDevices::TealGrn", 23)[6:22]
          )

# Get all unique nodes
all_nodes <- unique(c(data$source, data$target))

# Create graph
graph <- graph_from_data_frame(data, directed = FALSE, vertices = all_nodes)

# Add edge weights (probabilities)
E(graph)$weight <- data$prob

# Create layout
layout <- create_layout(graph, layout = "circle")

# Calculate label positions (slightly outside the circle)
label_radius <- 1.1
layout$label_x <- layout$x * label_radius
layout$label_y <- layout$y * label_radius
layout$label_angle <- atan2(layout$y, layout$x) * 180 / pi + 180

# Create plot
ggraph(layout) +
  geom_edge_arc(aes(edge_width = weight, alpha = weight), edge_colour = "gray50", strength = 0.2) +
  scale_edge_width(range = c(0, 2)) +
  scale_edge_alpha(range = c(0, 1)) +
  geom_node_point(size = 5, color = col2) +
  scale_color_manual(values = col2) +
  geom_node_text(aes(x = label_x, y = label_y, label = name, angle = label_angle), 
                 size = 4, hjust = 0.9, vjust = 0.9) +
  coord_fixed() +
  theme_void() +
  labs(title = "IFNG signaling - AD")



################

### Fig 2C

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1)]
order1 <- c("Harvard_C","Harvard_AD","Sichuan_C","Sichuan_AD","Chicago_C","Chicago_AD")
prot.combined$cohort.disease <- factor(prot.combined$cohort.disease, levels = order1)

my_comparisons <- list(c("Harvard_AD","Harvard_C"),c("Sichuan_AD","Sichuan_C"),
                       c("Chicago_AD","Chicago_C"))
Idents(prot.combined) <- "celltypes"
# Idents(prot.combined) <- "seurat_clusters"
gene1 <- "IL15"

prot.combined2 <- subset(prot.combined, subset = !(patient %in% c("SAMN33579879")))

p1 <- VlnPlot(prot.combined2, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors,
              idents = "Classical monocytes"  # "52"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.01, max(p1[[1]][["data"]][[gene1]])+1.5))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1
# ggsave(filename = paste0("vln_stats_IL15_mean_classic.mono_by.sex.age2.disease.pdf"),
#        plot = p1, width = 12, height = 10,
#        #dpi = 600,
#        device = "pdf")
# rm(p1)

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1)]
order1 <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = order1)

my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))
Idents(prot.combined) <- "seurat_clusters"
p1 <- VlnPlot(prot.combined, features = "IFNG",
              pt.size = 0.05, raster = F, group.by = "age.range.disease2", cols = colors,
              idents = "52"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IFNG"]])))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1


### Fig 3C

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1)]
order1 <- c("Harvard_C","Harvard_AD","Sichuan_C","Sichuan_AD","Chicago_C","Chicago_AD")
prot.combined$cohort.disease <- factor(prot.combined$cohort.disease, levels = order1)

my_comparisons <- list(c("Harvard_AD","Harvard_C"),c("Sichuan_AD","Sichuan_C"),
                       c("Chicago_AD","Chicago_C"))
Idents(prot.combined) <- "celltypes"
Idents(prot.combined) <- "seurat_clusters"
gene1 <- "IFNG"

prot.combined2 <- subset(prot.combined, subset = !(patient %in% c("SAMN33579879")))

p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors,
              idents =  "52"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.00000, max(p1[[1]][["data"]][[gene1]])+0.))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1



# by sex
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1)]
order1 <- c("male_C","male_AD","female_C","female_AD")
prot.combined$sex.disease <- factor(prot.combined$sex.disease, levels = order1)

my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))
Idents(prot.combined) <- "seurat_clusters"
p1 <- VlnPlot(prot.combined, features = "IFNG",
              pt.size = 0.05, raster = F, group.by = "sex.disease", cols = colors,
              idents = "62"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IFNG"]])))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1


# by APOE
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,2,2,2,1,1,1,2,1)]
order1 <- c("C_E2/E3","C_E3/E3","C_E3/E4","C_E4/E4","AD_E3/E3","AD_E3/E4","AD_E4/E4","C_","AD_")
prot.combined$geno.disease2 <- factor(prot.combined$geno.disease2, levels = order1)

my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))
Idents(prot.combined) <- "seurat_clusters"
p1 <- VlnPlot(prot.combined, features = "IFNG",
              pt.size = 0.05, raster = F, group.by = "geno.disease2", cols = colors,
              idents = "52"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IFNG"]])))#+
  # stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1

##### Fig 2D

# res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_NKT-like CD8+_DEGs.xlsx")
res <- read.xlsx("DEGs/age2.sex.cohort.regress/clusters/MAST_clusters_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_clus52_DEGs.xlsx")
colnames(res)[1] <- "genes"
hav <- c("STAT3","HLA-A","HLA-E","BTN3A2","CD3E","HLA-DRB1","HLA-DPB1","HLA-DPA1","SUMO1","MT2A","B2M","IFITM2","PSMB8")
EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-1.1, 1.1),
                ylim = c(0, 46),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "NKT-like CD8+ cluster3 - AD x C",
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




res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_NKT-like CD8+_DEGs.xlsx")
colnames(res)[1] <- "genes"
hav <- c("GZMB","NKG7","GZMA","SRGN","PRF1","B2M","CD2","CTSC","HCST","TUBB","FCGR3A")
EnhancedVolcano(res,
                lab = res$genes, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-1.1, 1),
                #ylim = c(0, 46),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "NKT-like CD8+ AD x C",
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


##### Fig 2E

file1 <- read.delim("GOBP_upregulated_NKTCD8_clus52_sel_reduced.csv", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3,4,5)
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
  ggtitle("NKT-like CD8+ cluster3 - AD x C") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()

file1 <- read.delim("GOBP_upregulated_NKTCD8_sel.csv", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3.6,4,5.5)
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
  ggtitle("NKT-like CD8+  AD x C") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()

##### Fig 2F and G
setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/")

file1 <- read.xlsx("all_mono_2_NK.NKTCD8_nichenet_vis_ligand_target4.xlsx") # for F
file1 <- read.xlsx("all_mono_2_effect.CD8_nichenet_vis_ligand_target4.xlsx") # for G
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
  ggtitle("Classical monocytes -> CD8+ T cells") +
  theme(plot.title = element_text(hjust = 0.3, size = 10)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$mean.interaction)), 
                       name = expression("Regulatory potential"#"Log"[2]*"(Fold-Enrichment)"
                       )) 

####### monocytes common signature brain and periphery
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

FeaturePlot(prot.combined, features = "Mono.markers1", raster = T, pt.size = 1.5,
            cols = c("gray80","#9E1021"),
            # split.by = "region_disease",
            # ncol = 2,
            min.cutoff = 0.5,
            max.cutoff = 2.,
            reduction = "umap") + theme(legend.position = "right")


hav <- read.xlsx("MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_Classical monocytes_DEGs.xlsx")
colnames(hav)[1] <- "genes"
hav2 <- read.xlsx("wilcox_myeloids_ADxHC_Monocytes . Macrophages_DEGs.xlsx")
colnames(hav2)[1] <- "genes"
hav3 <- hav2[hav2$genes %in% hav$genes,]
hav4 <- hav[hav$genes %in% hav3$genes,]

hav3 <- hav3[order(hav3$genes),]
hav4 <- hav4[order(hav4$genes),]

mydf1 <- data.frame(row.names = hav3$genes,
                   PBMC = hav4$avg_log2FC,
                   Brain = hav3$avg_log2FC)
                   
mydf1$group <- "Other"
mydf1$group[mydf1$PBMC > 0 & mydf1$Brain > 0] <- "Both_Positive"
mydf1$group[mydf1$PBMC < 0 & mydf1$Brain < 0] <- "Both_Negative"

selected_points <- c("IFNGR2","JAK1","IRF2","IRF3","PIAS1","PPARA","PPARG","CD36","TLR2","MYD88","TLR6","SREBF1","SREBF2",
                            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K4","MAP3K7","MAP3K14","MAP4K3","MAP4K5","MAPK8",
                            "CREB1","CREBBP",#,"CD86","LRP1","IRF1","IFNGR1","STAT1","MAP3K20","JAK2","MAP3K1","MAP3K3","MAP3K5","MAP4K4","MAPK1","TLR4",
                            "IL15","IL15RA")
label_data <- mydf1[rownames(mydf1) %in% selected_points, ]

ggscatter(mydf1, x = "PBMC", y = "Brain", #label = rownames(mydf1),
          color = "group", palette = c("Both_Positive" = "salmon", "Both_Negative" = "lightblue", "Other" = "grey70"),
          repel = T, #font.label = c(14, "bold.italic", "black"),
          add = "reg.line", # Add regression line
          add.params = list(color = "#AE123A", fill = "gray24", size = 1.4), # Customize reg. line
          conf.int = T ) + # Add confidence interval
          geom_text_repel(data = label_data, max.overlaps = 200,
                  aes(x = PBMC, y = Brain, label = rownames(label_data)),
                  fontface = "bold", size = 3) +
  # Add dashed lines at x=0 and y=0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
  stat_cor(method = "pearson", label.x = -1.3, label.y = 6) +
  annotate("text", x = -1, y = 5, label = (paste0("slope==", coef(lm(mydf1$Brain~mydf1$PBMC))[2])), parse = T) + 
  ggtitle("PBMC and Brain Monocytes correlation") + ylab(expression("AD x C Log"[2]*"(Fold-Change) - Brain")) + 
  xlab(expression("AD x C Log"[2]*"(Fold-Change) - PBMC"))

###############################        fig 3     ##############################
######## Fig 3A

res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_Classical monocytes_DEGs.xlsx")
hav <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","IRF2",#"IRF3","IRF1",
         "PIAS1","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
         #"MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
         "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
         "CREB1","CREBBP",#,"CD86"
         "IL15","IL15RA")

EnhancedVolcano(res,
                lab = res$gene, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-2, 2),
                # ylim = c(0, 46),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 6.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Classical monocytes - AD x C",
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

res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_Nonclassical monocytes_DEGs.xlsx")
hav <- c(#"IFNGR1",
         "IFNGR2","JAK1","JAK2","STAT1","IRF2",#"IRF3","IRF1",
         "PIAS1","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
         #"MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
         "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
         "CREB1","CREBBP",#,"CD86"
         "IL15","IL15RA")

EnhancedVolcano(res,
                lab = res$gene, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-2.5, 2),
                # ylim = c(0, 46),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 6.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Nonclassical monocytes - AD x C",
                subtitle = bquote(italic("")),
                pointSize = 2,
                shadeAlpha = 2,
                #lengthConnectors = unit(0.01, "npc"),
                #arrowheads = T,
                max.overlaps = 2000,
                #maxoverlapsConnectors = NULL,
                #min.segment.length = 0.0000001,
                #directionConnectors = "x",
                #parseLabels = FALSE,
                raster = FALSE,
                #typeConnectors = "open",
                #endsConnectors = "first",
                caption = "" #paste0("total = ", nrow(toptable), " variables"),
)

res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_Intermediate monocytes_DEGs.xlsx")
hav <- c("IFNGR1",#"IFNGR2","JAK2",
         "JAK1","IRF2",#"IRF3","IRF1","STAT1",
         "PIAS1","PPARG","CD36","TLR4",#"MYD88",#"TLR6","LRP1","SREBF1","SREBF2","TLR2","PPARA",
         #"MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", "MAPK1","MAP3K14","MAP4K4",
         "MAP3K7","MAP3K20","MAP4K3","MAP4K5","MAPK8",
         "CREB1",#"CREBBP",#,"CD86","IL15RA"
         "IL15")

EnhancedVolcano(zk.response0,
                lab = row.names(zk.response0), #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-2, 2),
                ylim = c(0, 75),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 6.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Intermediate monocytes - AD x C",
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

###### Fig 3B
#####  GO bubble

file1 <- read.delim("reactome_classical.mono_upregulated_sel_reduced.csv", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3.5,7,10.5)
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
  ggtitle("Intermediate monocytes AD upregulated genes") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       )) #+ 
#theme(panel.background = element_rect(fill = "white"))
#theme_classic() + 
#theme_minimal()


##### Fig 2E Dotplot

my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()


##### Fig 3C

iNK.react <- read.delim("iNK_up_AD_reactome.txt", header = T, skip = 11)
mNK.react <- read.delim("mNK_up_AD_reactome.txt", header = T, skip = 11)
NKT.CD8.react <- read.delim("NKT.CD8_up_AD_reactome.txt", header = T, skip = 11)

iNK.react2 <- iNK.react[iNK.react$upload_1..733. >= 9 & iNK.react$upload_1..fold.Enrichment. > 3 
                        #& iNK.react$upload_1..FDR. < 2.00e-02 
                          & iNK.react$Homo.sapiens...REFLIST..20580. <= 100 
                        #& iNK.react$upload_1..raw.P.value. < 2.00e-03
                        , ]

mNK.react2 <- mNK.react[mNK.react$upload_1..600. > 8 & mNK.react$upload_1..fold.Enrichment. > 3 
                        #& mNK.react$upload_1..FDR. < 2.00e-02 
                          & mNK.react$Homo.sapiens...REFLIST..20580. <= 100 
                        #& mNK.react$upload_1..raw.P.value. < 2.00e-03
                        , ]

NKT.CD8.react2 <- NKT.CD8.react[NKT.CD8.react$upload_1..610. >= 9 & NKT.CD8.react$upload_1..fold.Enrichment. > 3 
                         #       & NKT.CD8.react$upload_1..FDR. < 2.00e-02 
                                &  NKT.CD8.react$Homo.sapiens...REFLIST..20580. <= 100 
                          #      & NKT.CD8.react$upload_1..raw.P.value. < 2.00e-03
                                , ]

common_elements <- Reduce(intersect, list(iNK.react2$Reactome.pathways, mNK.react2$Reactome.pathways, NKT.CD8.react2$Reactome.pathways))

iNK.react3 <- iNK.react[iNK.react$Reactome.pathways %in% common_elements,]
mNK.react3 <- mNK.react[mNK.react$Reactome.pathways %in% common_elements,]
NKT.CD8.react3 <- NKT.CD8.react[NKT.CD8.react$Reactome.pathways %in% common_elements,]

write.xlsx(as.data.frame(NKT.CD8.react3), rowNames = T, file="NKT.CD8.up.AD.reactome3.xlsx")

sel1 <- read.delim("NKT.CD8_up_reactome_sel2.txt", header = T)

iNK.react4 <- iNK.react[iNK.react$Reactome.pathways %in% sel1$Reactome.pathways,]
mNK.react4 <- mNK.react[mNK.react$Reactome.pathways %in% sel1$Reactome.pathways,]
NKT.CD8.react4 <- NKT.CD8.react[NKT.CD8.react$Reactome.pathways %in% sel1$Reactome.pathways,]

combined_df <- bind_rows(
  iNK.react4 %>% mutate(source = "iNK"),
  mNK.react4 %>% mutate(source = "mNK"),
  NKT.CD8.react4 %>% mutate(source = "NKT-like CD8+")
) %>%
  mutate(log10_value1 = -log10(upload_1..raw.P.value.)) %>% mutate(enrichment = as.integer(upload_1..fold.Enrichment.))

coluna1 <- sel1$Reactome.pathways[15:1]
mid <- 0

# Create the plot
ggplot() +
  geom_point(data=combined_df, aes(x = source, y = factor(Reactome.pathways, levels = coluna1),
                                   size = log10_value1, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3,4.5,7)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    axis.text.x = element_text(face="bold", color = "black",size=12, angle = 45, hjust = 0.9), 
    axis.text.y = element_text(face="bold",color = "black",size=12),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Common enriched pathways of upregulated DEGs in iNK, mNK and NKT-like CD8+") +
  theme(plot.title = element_text(hjust = 0.9)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(combined_df$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       ))















################

### Fig 2F

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

my_comparisons <- c("male_C","male_AD","female_C","female_AD")
prot.combined$sex.disease <- factor(prot.combined$sex.disease, levels = my_comparisons)

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "sex.disease",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

### Fig suppl 2
Idents(prot.combined) <- "celltypes"

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

# my_comparisons <- c("male_C","male_AD","female_C","female_AD")
# prot.combined$sex.disease <- factor(prot.combined$sex.disease, levels = my_comparisons)
prot.harvard <- subset(prot.combined, subset = cohort %in% c("Harvard"))
prot.sichuan <- subset(prot.combined, subset = cohort %in% c("Sichuan"))
prot.chicago <- subset(prot.combined, subset = cohort %in% c("Chicago"))

DotPlot(prot.chicago, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "cohort.disease",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()



###########

##### Fig 2G

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1,2,1)]
prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")
order1 <- c("50-70_male_C","50-70_male_AD",
            "50-70_female_C","50-70_female_AD",
            "71-90_male_C","71-90_male_AD",
            "71-90_female_C","71-90_female_AD")
prot.combined$sex.age2.disease <- factor(prot.combined$sex.age2.disease, levels = order1)

my_comparisons <- list(c("50-70_male_AD","50-70_male_C"),c("71-90_male_AD","71-90_male_C"),
                       c("50-70_female_AD","50-70_female_C"),c("71-90_female_AD","71-90_female_C"))

p1 <- VlnPlot(prot.combined, features = "TSPO",
              pt.size = 0.05, raster = F, group.by = "sex.age2.disease", cols = colors,
              idents = "Intermediate monocytes"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["TSPO"]])+1.8))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1
ggsave(filename = paste0("vln_stats_IL15_mean_classic.mono_by.sex.age2.disease.pdf"),
       plot = p1, width = 12, height = 10,
       #dpi = 600,
       device = "pdf")
rm(p1)

# by APOE
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,2,2,2,1,1,1,2,1)]
order1 <- c("C_E2/E3","C_E3/E3","C_E3/E4","C_E4/E4","AD_E3/E3","AD_E3/E4","AD_E4/E4","C_","AD_")
prot.combined$geno.disease2 <- factor(prot.combined$geno.disease2, levels = order1)

# my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))
Idents(prot.combined) <- "celltypes"
p1 <- VlnPlot(prot.combined, features = "IL15",
              pt.size = 0.05, raster = F, group.by = "geno.disease2", cols = colors,
              idents = "Classical monocytes"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][["IL15"]])))#+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1




###############################        fig 4     ##############################
###### Fig 3A

res <- read.xlsx("DEGs/age2.sex.cohort.regress/celltypes/MAST_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_Effector CD8+ T cells_DEGs.xlsx")
hav <- c("TIGIT","LTA","LTB","SLAMF6","IL2RB","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","IL2RG",
         "PRF1","GZMK","GZMH","GZMA", #"GZMB","GZMM",
         "NKG7","TYROBP",#"FASLG","KLRB1",
         "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
         "ENTPD1")

EnhancedVolcano(res,
                lab = res$gene, #NA,
                selectLab = hav,#tmarkers[1,],
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.1,
                FCcutoff = 0.05,
                xlim = c(-1, 0.75),
                ylim = c(0, 105),
                drawConnectors = T,
                #widthConnectors = 0.8,
                labSize = 5.0,
                labFace = 'bold',
                #boxedLabels = T,
                title = "Effector CD8+ T cells - AD x C",
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



###### Fig 3A
#####  GO bubble

file1 <- read.delim("GO_BP_effector_t.cd8_upregulated_sel2.csv", header = T)
colunas1 <- file1$Term[23:1]
mid <- 0

ggplot() + geom_point(data=file1, aes(x = "",#genes, 
                                      y = factor(Term, levels = colunas1), 
                                      size = log10pval, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(3.5,5,6.5)
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
  ggtitle("Effector CD8+ T cells upregulated genes") +
  theme(plot.title = element_text(hjust = 1)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(file1$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       ))

##### Fig 3B Dotplot

my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)

genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
            "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
            "NKG7","TYROBP",#"FASLG","KLRB1",
            "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
            "ENTPD1"
            #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
            
) # effector cd8
Idents(prot.combined) <- "celltypes"
DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()


####
VlnPlot(prot.combined, )

### Fig 3C
my_comparisons <- c("male_C","male_AD","female_C","female_AD")
prot.combined$sex.disease <- factor(prot.combined$sex.disease, levels = my_comparisons)
genes1 <- c(#"CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
  "SLAMF6","IRF2","GNLY","TOX","IFNG","STAT4","NKG7","GZMH","GZMA", #males
  "GZMK","PRF1","TIGIT","LTA","IL2RB","IL2RG","EOMES","TYROBP","BATF","CD244","ENTPD1"#females
) # effector cd8
genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
  "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
  "NKG7","TYROBP",#"FASLG","KLRB1",
  "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
  "ENTPD1"
  #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
  
) # effector cd8
DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "sex.disease",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()



################              APOE genotype




prot.combined2 <- subset(prot.combined, subset = !(geno.disease2 %in% c("AD_","C_","C_E2/E3" ,"C_E3/E4","AD_E3/E4" #,"C_E3/E3","C_E3/E4","C_E4/E4"
                                                                        )) )
my_comparisons <- c("AD_E3/E3","AD_E3/E4","AD_E4/E4")
my_comparisons <- c("C_E3/E3","C_E4/E4","AD_E3/E3","AD_E4/E4")
prot.combined2$geno.disease2 <- factor(prot.combined2$geno.disease2, levels = my_comparisons)
prot.combined2$age.geno <- paste(prot.combined2$age.range2, prot.combined2$geno, sep = "_")
prot.combined2$sex.geno <- paste(prot.combined2$sex, prot.combined2$geno, sep = "_")

IL15.sig <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
              "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
              "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
              "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
              "CREB1","CREBBP",
              "IL15","IL15RA"#,"CD86"
) # classic Monocytes

prot.combined2 <- AddModuleScore(
  prot.combined2,
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
colors <- brewer.pal(n=6,name = "Dark2")
p1 <- VlnPlot(prot.combined2, features = gene1,
              pt.size = 0.05, raster = F, group.by = "age.geno", cols = colors,
              idents = "Classical monocytes"
              # "CD8+ T cells"
              # "Homeostatic Microglia"
              # "DAM Microglia"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) #+
# stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1






Idents(prot.combined2) <- "seurat_clusters"
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(1,1,1)]
gene1 <- "IFNG"
my_comparisons <- list(c("AD_E3/E3","AD_E3/E4"),c("AD_E3/E3","AD_E4/E4"),c("AD_E3/E4","AD_E4/E4"))
my_comparisons <- list(c("E3/E3","E3/E4"),c("E3/E3","E4/E4"),c("E3/E4","E4/E4"))
p1 <- VlnPlot(prot.combined2, features = gene1, pt.size = 0.01, raster = F,
              # ncol = 2, stack = T, flip = T
              group.by = "geno", cols = colors,
              idents = c("52") #,"41,"52","62"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[gene1]])+1.5)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1

Idents(prot.combined2) <- "celltypes"
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(1,1,1,1,1,1)]
gene1 <- "IL15"
my_comparisons <- list(c("male_E3/E3","male_E3/E4"),c("male_E3/E3","male_E4/E4"),c("male_E3/E4","male_E4/E4"),
                       c("female_E3/E3","female_E3/E4"),c("female_E3/E3","female_E4/E4"),c("female_E3/E4","female_E4/E4"))
p1 <- VlnPlot(prot.combined2, features = gene1, pt.size = 0.01, raster = F,
              # ncol = 2, stack = T, flip = T
              group.by = "sex.geno", cols = colors,
              idents = c("Classical monocytes") #,"41,"52","62"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+2.5)) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.05
p1

genes1 <- c(#"CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
  "SLAMF6","IRF2","GNLY","TOX","IFNG","STAT4","NKG7","GZMH","GZMA", #males
  "GZMK","PRF1","TIGIT","LTA","IL2RB","IL2RG","EOMES","TYROBP","BATF","CD244","ENTPD1"#females
) # effector cd8
genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
  "PRF1","GZMH","GZMA", "GZMB","GZMM",#"GZMK",
  "NKG7","TYROBP",#"FASLG","KLRB1",
  "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
  "ENTPD1"
  #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
) # effector cd8

genes1 <- c("CD8B","CD8A","IRF2","IRF3",  ## for genotype  "IL2RA","TNF","IRF1","IRF4","GZMB","KLRB1","TBX21","FYN","CD28","NFATC1","NR4A1","STAT4","IFNG",
            "GZMM","FASLG","ZAP70","LCK","TCF7","NFATC2","LTA",
  "SLAMF6","BATF","CD244","IL2RG","GNLY","TOX","NKG7","GZMH","GZMA", #
  "PRF1","TIGIT","IL2RB","EOMES","TYROBP" #,"ENTPD1"
) # effector cd8

DotPlot(prot.combined2, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 10,
        # cluster.idents = T,
        #scale = F,
        group.by = "geno.disease2",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()

genes1 <- c("IFNGR1","IFNGR2","JAK1","STAT2","STAT3","IRF1", # "IRF2","IRF3","PIAS1","RUNX3","CD36","MYD88","TLR6","SREBF1",  ### for genotype
            "PPARA","PPARG","TLR2","TLR4","LRP1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5",  
            "MAP3K7","MAP3K20","MAP4K4","MAPK8",  # "MAP3K14","MAP4K3","MAP4K5","MAPK1",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86","JAK2","STAT1",
) # classic Monocytes

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

DotPlot(prot.combined2, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "geno.disease2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

# prot.combined2M <- subset(prot.combined2, subset = sex %in% c("male"))
# prot.combined2F <- subset(prot.combined2, subset = sex %in% c("female"))
# prot.combined2Late <- subset(prot.combined2, subset = age.range2 %in% c("71-90"))
# prot.combined2Early <- subset(prot.combined2, subset = age.range2 %in% c("50-70"))
prot.combined2C <- subset(prot.combined2, subset = disease %in% c("C"))
prot.combined2AD <- subset(prot.combined2, subset = disease %in% c("AD"))

# my_comparisons <- c("C_E2/E3","C_E3/E3","C_E3/E4","C_E4/E4",
#                     "AD_E3/E3","AD_E3/E4","AD_E4/E4")
# prot.combined2F$geno.disease2 <- factor(prot.combined2F$geno.disease2, levels = my_comparisons)

genes1 <- c(#"CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
  "SLAMF6","IRF2","GNLY","TOX","IFNG","STAT4","NKG7","GZMH","GZMA", #males
  "GZMK","PRF1","TIGIT","LTA","IL2RB","IL2RG","EOMES","TYROBP","BATF","CD244","ENTPD1"#females
) # effector cd8
genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
  "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
  "NKG7","TYROBP",#"FASLG","KLRB1",
  "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
  "ENTPD1"
  #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
) # effector cd8

genes1 <- c("CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1",
               "SLAMF6","IRF2","GNLY","TOX","IFNG","STAT4","NKG7","GZMH","GZMA", #males
               "PRF1","TIGIT","LTA","IL2RB","IL2RG","EOMES","TYROBP","BATF","CD244","ENTPD1"#females
             ) # effector cd8

DotPlot(prot.combined2, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 10,
        # cluster.idents = T,
        #scale = F,
        group.by = "geno.disease2",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

DotPlot(prot.combined2Late, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 10,
        # cluster.idents = T,
        #scale = F,
        group.by = "geno.disease2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()


VlnPlot(prot.combined2Early, features = "IL15", idents = "Classical monocytes",
        group.by = "geno.disease2")

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,2,2,2,1,1,1)]

p1 <- VlnPlot(prot.combined2F, features = "IL15",
              pt.size = 0.05, raster = F, group.by = "geno.disease2", cols = colors,
              idents = c("Classical monocytes")
) + theme(legend.position = "none")+ stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95)
p1$layers[[2]]$aes_params$alpha <- 0.1
p1

p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.00001, max(p1[[1]][["data"]][["IL15"]])+0))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1


# prot.combined3 <- subset(prot.combined, subset = cohort2 %in% c("BWH"))
# prot.combined3$mmse <- prot.combined$age
# prot.combined3$mmse.range <- prot.combined$age
# age1 <- c("68","74","75","78","81")
# mmse.range <- c("Mild","Normal","Normal","Mild","Normal")
# mmse <- c("20","29","25","24","27")
# 
# for (i in 1:length(age1)) {
#   prot.combined3$mmse <- recode(prot.combined3$mmse, "age1[i] = mmse[i]")
#   prot.combined3$mmse.range <- recode(prot.combined3$mmse.range, "age1[i] = mmse.range[i]")
# }
# 
# genes1 <- c(#"CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
#   "SLAMF6","IRF2","GNLY","TOX","IFNG","STAT4","NKG7","GZMH","GZMA", #males
#   "GZMK","PRF1","TIGIT","LTA","IL2RB","IL2RG","EOMES","TYROBP","BATF","CD244","ENTPD1"#females
# ) # effector cd8
# genes1 <- c( #"TNF","CD8B","CD8A","IL2RA",
#   "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES","IRF2","GNLY", #"IRF1","IRF3","IRF4","LTB"
#   "PRF1","GZMK","GZMH","GZMA", "GZMB",#"GZMM",
#   "NKG7","TYROBP",#"FASLG","KLRB1",
#   "BATF","TOX","IFNG","STAT4","CD244",#"TBX21",
#   "ENTPD1"
#   #"IL2RB","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1"
#   
# ) # effector cd8
# 
# DotPlot(prot.combined3, features = rev(genes1),
#         cols = "RdBu",
#         col.max = 20,
#         dot.scale = 10,
#         # cluster.idents = T,
#         #scale = F,
#         group.by = "mmse.range",
#         idents = "Effector CD8+ T cells",
# ) + RotatedAxis() + coord_flip()
# 
# genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
#             "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
#             "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
#             "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
#             "CREB1","CREBBP",
#             "IL15","IL15RA"#,"CD86"
# ) # classic Monocytes
# 
# DotPlot(prot.combined3, features = rev(genes1),
#         cols = "RdBu",
#         col.max = 20,
#         dot.scale = 10,
#         # cluster.idents = T,
#         #scale = F,
#         group.by = "mmse",
#         idents = "Classical monocytes",
# ) + RotatedAxis() + coord_flip()
# 





#### Fig 3X

celltypes <- c(
  "Platelets_C",
  "Platelets_AD",
  "Classical monocytes_C",
  "Classical monocytes_AD",
  "Intermediate monocytes_C",
  "Intermediate monocytes_AD",
  "Nonclassical monocytes_C",
  "Nonclassical monocytes_AD",
  "mo-DC_C",
  "mo-DC_AD",
  "pDC_C",
  "pDC_AD",
  "iNK_C",
  "iNK_AD",
  "mNK_C",
  "mNK_AD",
  "NKT-like CD4+_C",
  "NKT-like CD4+_AD",
  "NKT-like CD8+_C",
  "NKT-like CD8+_AD",
  "Memory T CD8+_C",
  "Memory T CD8+_AD",
  "Naive CD8+ T cells_C",
  "Naive CD8+ T cells_AD",
  "Effector CD8+ T cells_C",
  "Effector CD8+ T cells_AD",
  "Exhausted CD8+ T cells_C",
  "Exhausted CD8+ T cells_AD",
  "Naive CD4+ T cells_C",
  "Naive CD4+ T cells_AD",
  "Memory T CD4+_C",
  "Memory T CD4+_AD",
  "Th1_C",
  "Th1_AD",
  "Cytotoxic CD4+ T cells_C",
  "Cytotoxic CD4+ T cells_AD",
  "Treg_C",
  "Treg_AD",
  "gamma-delta T cells_C",
  "gamma-delta T cells_AD",
  "Naive B cells_C",
  "Naive B cells_AD",
  "Memory B cells_C",
  "Memory B cells_AD",
  "Cytotoxic B cells_C",
  "Cytotoxic B cells_AD",
  "mature B cells_C",
  "mature B cells_AD",
  "Plasma B cells_C",
  "Plasma B cells_AD"
)

prot.combined$celltypes.disease <- factor(prot.combined$celltypes.disease, levels = celltypes)

# prot.combined$disease <- factor(prot.combined$disease, levels = c("C","AD"))
Idents(prot.combined) <- "celltypes.disease"
DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        #split.by = "disease",
        idents = c("Naive CD8+ T cells_C","Naive CD8+ T cells_AD",
                   "Memory T CD8+_C","Memory T CD8+_AD",
                   #"Exhausted CD8+ T cells_C","Exhausted CD8+ T cells_AD",
                   "Effector CD8+ T cells_C","Effector CD8+ T cells_AD"),
) + RotatedAxis() + coord_flip()

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1)]

my_comparisons <- list(c("Naive CD8+ T cells_C","Naive CD8+ T cells_AD"),
                       c("Memory T CD8+_C","Memory T CD8+_AD"),
                       c("Effector CD8+ T cells_C","Effector CD8+ T cells_AD"))

p1 <- VlnPlot(prot.combined, features = "IL2RG",
              pt.size = 0.05, raster = F, group.by = "celltypes.disease", cols = colors,
              idents = c("Naive CD8+ T cells_C","Naive CD8+ T cells_AD",
                          "Memory T CD8+_C","Memory T CD8+_AD",
                          #"Exhausted CD8+ T cells_C","Exhausted CD8+ T cells_AD",
                          "Effector CD8+ T cells_C","Effector CD8+ T cells_AD")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 30, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.00001, max(p1[[1]][["data"]][["IL2RG"]])+0))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.1
p1


############ Fig 3L

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/public_data/AD/CSF/Seurat/")

prot.combined2 <- readRDS("./obj_unintegrated.Rds")
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

DimPlot(prot.combined2, reduction = "umap.unintegrated", 
        label = T, raster=F, repel = T, group.by = "seurat_clusters"
)

########## 
iNK.react <- read.delim("GOBP_up_ADxC_effector_t.cd8.txt", header = T, skip = 11)
mNK.react <- read.delim("GOBP_up_ADxC_memory_cd8.txt", header = T, skip = 11)
NKT.CD8.react <- read.delim("GOBP_up_ADxC_naive_cd8.txt", header = T, skip = 11)

iNK.react2 <- iNK.react[iNK.react$upload_1..206. >= 9 & iNK.react$upload_1..fold.Enrichment. > 3 
                        #& iNK.react$upload_1..FDR. < 2.00e-02 
                        #& iNK.react$Homo.sapiens...REFLIST..20580. <= 100 
                        #& iNK.react$upload_1..raw.P.value. < 2.00e-03
                        , ]

mNK.react2 <- mNK.react[mNK.react$upload_1..406. > 9 & mNK.react$upload_1..fold.Enrichment. > 3 
                        #& mNK.react$upload_1..FDR. < 2.00e-02 
                        #& mNK.react$Homo.sapiens...REFLIST..20580. <= 100 
                        #& mNK.react$upload_1..raw.P.value. < 2.00e-03
                        , ]

NKT.CD8.react2 <- NKT.CD8.react[NKT.CD8.react$upload_1..445. >= 9 & NKT.CD8.react$upload_1..fold.Enrichment. > 3 
                                #       & NKT.CD8.react$upload_1..FDR. < 2.00e-02 
                                #&  NKT.CD8.react$Homo.sapiens...REFLIST..20580. <= 100 
                                #      & NKT.CD8.react$upload_1..raw.P.value. < 2.00e-03
                                , ]
# iNK.react2 <- iNK.react
# mNK.react2 <- mNK.react
# NKT.CD8.react2 <- NKT.CD8.react

common_elements <- Reduce(intersect, list(iNK.react2$GO.biological.process.complete, mNK.react2$GO.biological.process.complete, NKT.CD8.react2$GO.biological.process.complete))

iNK.react3 <- iNK.react[iNK.react$GO.biological.process.complete %in% common_elements,]
mNK.react3 <- mNK.react[mNK.react$GO.biological.process.complete %in% common_elements,]
NKT.CD8.react3 <- NKT.CD8.react[NKT.CD8.react$GO.biological.process.complete %in% common_elements,]

write.xlsx(as.data.frame(NKT.CD8.react3), rowNames = T, file="Naive.CD8.up.AD.GOBP3.xlsx")

sel1 <- read.delim("Naive.CD8.up.AD.GOBP_sel.txt", header = T)

iNK.react4 <- iNK.react[iNK.react$GO.biological.process.complete %in% sel1$GO.biological.process.complete,]
mNK.react4 <- mNK.react[mNK.react$GO.biological.process.complete %in% sel1$GO.biological.process.complete,]
NKT.CD8.react4 <- NKT.CD8.react[NKT.CD8.react$GO.biological.process.complete %in% sel1$GO.biological.process.complete,]

iNK.react4$upload_1..fold.Enrichment. <- as.numeric(as.character(iNK.react4$upload_1..fold.Enrichment.))
NKT.CD8.react4$upload_1..fold.Enrichment. <- as.numeric(as.character(NKT.CD8.react4$upload_1..fold.Enrichment.))
mNK.react4$upload_1..fold.Enrichment. <- as.numeric(as.character(mNK.react4$upload_1..fold.Enrichment.))

combined_df <- bind_rows(
  iNK.react4 %>% mutate(source = "Effector T CD8+"),
  mNK.react4 %>% mutate(source = "Memory T  CD8+"),
  NKT.CD8.react4 %>% mutate(source = "Naive T CD8+")
) %>%
  mutate(log10_value1 = -log10(upload_1..raw.P.value.)) %>% mutate(enrichment = as.integer(upload_1..fold.Enrichment.))

coluna1 <- sel1$GO.biological.process.complete[15:1]
mid <- 0

# Create the plot
ggplot() +
  geom_point(data=combined_df, aes(x = source, y = factor(GO.biological.process.complete, levels = coluna1),
                                   size = log10_value1, fill = enrichment), alpha = 1, shape = 21) +
  scale_size(range = c(6, 12), name = expression("-Log"[10]*"("*italic("p")*"-value)"), breaks = c(4,5,7)
  ) + 
  xlab("") + #xlim(0,15) +
  ylab("") + theme(#axis.title.x = element_text(face="bold", color = "black", size=12),
    axis.text.x = element_text(face="bold", color = "black",size=12, angle = 45, hjust = 0.9), 
    axis.text.y = element_text(face="bold",color = "black",size=12),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "black", size = 0.1),
    legend.position = "right") +
  ggtitle("Common enriched pathways of upregulated DEGs in CD8+ T cells") +
  theme(plot.title = element_text(hjust = 0.9)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)* max(abs(combined_df$enrichment)), 
                       name = expression("Fold-Enrichment"#"Log"[2]*"(Fold-Enrichment)"
                       ))












######### Fig 3M
my_comparisons <- c("HC","MCI","AD")
prot.combined2$disease <- factor(prot.combined2$disease, levels = my_comparisons)

genes1 <- c(#"CD8B","CD8A","IL2RA","TNF","IRF1","IRF3","IRF4","GZMB","GZMM","IFNG","STAT4","CD244","ENTPD1","SLAMF6",
  #"FASLG","KLRB1","TBX21","ZAP70","LCK","FYN","TCF7","CD28","NFATC1","NFATC2","NR4A1","NKG7","TIGIT","LTA",
  "IRF2","BATF","GNLY","TOX","GZMH","GZMK","GZMA","PRF1","IL2RB","IL2RG","EOMES",
  "TYROBP","FASLG","PDCD1","LTB","IRF1","CD8A","CD8B"
) # effector cd8

DotPlot(prot.combined2, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "disease",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()

######### Sup Fig 3B
my_comparisons <- c("HC","MCI","AD")
prot.combined2$disease <- factor(prot.combined2$disease, levels = my_comparisons)

genes1 <- c("IFNGR1","IFNGR2",#"JAK1","JAK2","STAT1","STAT2","STAT3","IRF2","IRF3","IRF1","PIAS1","RUNX3","PPARG",
            "PPARA","TLR2","TLR4",#"TLR6","LRP1","SREBF1","SREBF2","MYD88","CD36",
            "MAPK14","MAP2K1","MAP2K5", #"MAP2K4","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5","MAP4K3","MAP4K5",
            "MAP3K7","MAP3K14","MAP3K20","MAP4K4","MAPK1","MAPK8",
            #"CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic mono

DotPlot(prot.combined2, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "disease",
        idents = "Classical monocytes"
) + RotatedAxis() + coord_flip()


##### Fig 3N

colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,3,1)]
my_comparisons <- c("HC","MCI","AD")
prot.combined2$disease <- factor(prot.combined2$disease, levels = my_comparisons)

my_comparisons <- list(c("AD","HC"),c("AD","MCI"),c("MCI","HC"))

p1 <- VlnPlot(prot.combined2, features = "IFNG",
              pt.size = 0.05, raster = F, group.by = "disease", cols = colors,
              #idents = "Effector CD8+ T cells"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.0001, max(p1[[1]][["data"]][["IFNG"]])+1.2))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.3
p1
ggsave(filename = paste0("vln_stats_IL15_mean_classic.mono_by.sex.age2.disease.pdf"),
       plot = p1, width = 12, height = 10,
       #dpi = 600,
       device = "pdf")
rm(p1)



###############################      OLeg P01

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/PPG_Oleg_AD/")

my_comparisons <- c("50-70_C","50-70_AD","71-90_C","71-90_AD")
prot.combined$age.range.disease2 <- factor(prot.combined$age.range.disease2, levels = my_comparisons)

genes1 <- c("HLA-DRB1","CD33","GRN","MS4A6A","FCER1G","BIN1","TREM1","CCR2","PSEN1","HAVCR2","MAPT","APP","GSAP","BACE1" #,"RORA","APOE","TREM2"
) # classic Monocytes

genes1 <- c("HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","FCER1G","CTSH","CTSB","CD33","GRN","ABI3","BCKDK","LILRB2","SPI1","PILRA","ZYX",   
            "NDUFS2","RIN3","MS4A6A","SCARB2","CCR2",
            "PSEN1","APP","HAVCR2","SCIMP","SNX1","PICALM"
            #"BIN1","TREM1","MAPT","GSAP","BACE1"#,"RORA","APOE","TREM2"
) 

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "disease",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

DotPlot(prot.combined, features = ROS,
        # cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()














DotPlot(myeloids, features = genes1,
        col.max = 20,
        cols = "RdBu",
        dot.scale = 10, 
        cluster.idents = F, group.by = "celltype.disease",
        idents = c(#"Microglia M0_C",
          #"Microglia M0_AD",
          #"Microglia MgND_C",
          #"Microglia MgND_AD",
          "Monocytes / Macrophages_C",
          "Monocytes / Macrophages_AD"),
        #scale = F,
        #split.by = "disease"
) + RotatedAxis() + coord_flip()


prot.combined2$celltype.disease <- paste(prot.combined2$celltypes, prot.combined2$disease, sep = "_")
Idents(prot.combined2) <- "celltype.disease"

celltypes2 <- c(
  "Microglia M0_C",
  "Microglia M0_AD",
  "Microglia MgND_C",
  "Microglia MgND_AD",
  "Monocytes / Macrophages_C",
  "Monocytes / Macrophages_AD",
  "CD8+ T cells_C",
  "CD8+ T cells_AD"
)

prot.combined2$celltype.disease <- factor(prot.combined2$celltype.disease, levels = celltypes2)


my_comparisons <- list(c("AD","C"))
Idents(myeloids) <- "celltypes"

for (i in 1:length(genes1)){
  p1 <- VlnPlot(myeloids, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "disease",
                idents = "Monocytes / Macrophages"
  ) + theme(legend.position = "none")
  p1 <- p1 + stat_summary(fun = median, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("Fig1_P01_vln_stats_",genes1[i],"_mean_mono.mac_brain_by.disease.png"),
         plot = p1, width = 3, height = 8,
         #dpi = 600,
         device = "png")
  rm(p1)
}

FeaturePlot(prot.combined, features = #c("HLA-DRB1","CD33","GRN"), 
              c(#"MS4A6A","FCER1G","BIN1"#,
                #"TREM1","CCR2","PSEN1"#,
                #"HAVCR2","MAPT","APP"#,
                #"GSAP","BACE1" #,"RORA","APOE","TREM2"
                "CCR2"
              ),
            pt.size = 2,
            #min.cutoff = 0.5,
            max.cutoff = 1.5,
            #ncol = 3, 
            reduction = "umap", raster = T, 
            split.by = "disease"
) + theme(legend.position = "right")



genes1 <- c("CCR2")
my_comparisons <- list(c("AD","C"))
for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "disease",
                idents = "Classical monocytes"
  ) + theme(legend.position = "none")
  p1 <- p1 + stat_summary(fun = median, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.0001, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.1
  ggsave(filename = paste0("P01_vln_stats_",genes1[i],"_mean_classic.mono_PBMC_by.disease.png"),
         plot = p1, width = 3, height = 8,
         #dpi = 600,
         device = "png")
  rm(p1)
}

my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "age.range.disease2",
                idents = "Classical monocytes"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("P01_vln_stats_",genes1[i],"_mean_classic.mono_by.age.disease2.pdf"),
         plot = p1, width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}














##################################################      Supplementary    ##########################


######### DEGs by cohort

setwd("/media/patrick/JAMBERT/Bioinfo/weiner_lab/izzy/scRNAseq/human/Seurat_integrated_AD2-3/BPCells2/DEGs/by.cohort/")

celltypes <- c(
  "Platelets",
  "Classical monocytes",
  "Intermediate monocytes",
  "Nonclassical monocytes",
  "mo-DC",
  "pDC",
  "iNK",
  "mNK",
  "NKT-like CD4+",
  "NKT-like CD8+",
  "Effector CD8+ T cells",
  "Naive CD8+ T cells",
  "Memory T CD8+",
  "Naive CD4+ T cells",
  "Memory T CD4+",
  "Th1",
  "Cytotoxic CD4+ T cells",
  "Treg",
  "gamma-delta T cells",
  "Naive B cells",
  "Memory B cells",
  "mature B cells",
  "Plasma B cells",
  "Cytotoxic B cells"#,
  # "Exhausted CD8+ T cells"
)

cohort <- c("Chicago","Harvard","Sichuan")

prot.combined$celltype.cohort.disease <- paste(prot.combined$celltypes, prot.combined$cohort.disease, sep = "_")
Idents(prot.combined) <- "celltype.cohort.disease"

for (i in 1:length(celltypes)) {
  for (e in 1:length(cohort)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_",cohort[e], "_AD"), 
                              ident.2 = paste0(celltypes[i],"_",cohort[e], "_C"),
                              slot = "data",
                              assay = "RNA",
                              features = NULL,
                              logfc.threshold = 0,
                              test.use = "wilcox",
                              min.pct = 0.3,
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
    zk.response0 <- zk.response0[zk.response0$p_val_adj < 0.1,]
    #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("wilcox_RNA_",cohort[e],"_ADxHC_", celltypes[i],"_DEGs.xlsx"))
    rm(zk.response0)
  }
}

for (i in 1:length(celltypes)) {
  for (e in 1:length(cohort)) {
    zk.response0 <- FindMarkers(prot.combined, ident.1 = paste0(celltypes[i],"_",cohort[e], "_AD"), 
                                ident.2 = paste0(celltypes[i],"_",cohort[e], "_C"),
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
                                latent.vars = c("age.range2","sex"),
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
    #write.csv(as.data.frame(zk.response0),file=paste0("wilcox_ADxHC_", cell.types[i], "_DEGs.csv"))
    write.xlsx(as.data.frame(zk.response0), rowNames = T,file=paste0("MAST_RNA_",cohort[e],"_ADxHC_", celltypes[i],"_DEGs.xlsx"))
    rm(zk.response0)
  }
}

# DEGs by cell type
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
  # "Exhausted CD8+ T cells",
)

data <- data.frame(matrix(ncol = 2,nrow = 0))

for (i in 1:length(celltypes)){
  list1 <- read.xlsx(paste0("MAST_RNA_Sichuan_ADxHC_",celltypes[i],"_DEGs.xlsx"))
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
new_row <- data.frame(clusters="Exhausted CD8+ T cells", DEGs=0,upregulated=0,downregulated=0)
new_df <- rbind(data[1:13, ], new_row, data[(14:24), ])
coluna1 <- new_df$clusters
coluna1 <- rev(coluna1) 
new_df$clusters <- factor(new_df$clusters, levels=coluna1)


axis_margin <- 5.5

p1 <- ggplot(new_df, aes(y=factor(clusters,levels=coluna1), x=downregulated)) +
  geom_col(fill = "slategray1") +
  geom_text(aes(label = downregulated), hjust = 1.5) + 
  scale_x_reverse(
    limits = c(700, 0),
    breaks = seq(500, 0, by = -100)
  ) + 
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
    axis.text.y = element_blank(),
  )

p2 <- ggplot(new_df, aes(y=factor(clusters,levels=coluna1), x=upregulated)) +
  geom_col(fill = "salmon") +
  geom_text(aes(label = upregulated), hjust = -0.5) + 
  scale_x_continuous(
    limits = c(0, 500),
    breaks = seq(0, 500, by = 100),
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
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
  )

ggarrange(p1, p2)



############ GWAS

hav <- t(unique(read.delim("AD-GWAS.txt", header = F)))
strips <- "^\\s+|\\s+$"
hav <- gsub(strips,"",hav)

Idents(prot.combined) <- "celltypes.disease"

avg.exp <- AverageExpression(
  prot.combined,
  assays = "RNA",
  features = hav,
  return.seurat = FALSE,
  group.by = "celltypes.disease",
  add.ident = NULL,
  layer = "data",
  # slot = deprecated(),
  verbose = TRUE
)

write.xlsx(as.data.frame(avg.exp$RNA), rowNames = T,file=paste0("AD-GWAS_genes_avg.exp_PBMC.xlsx"))

hav2 <- read.xlsx("AD-GWAS_genes_avg.exp_PBMC.xlsx")
hav2 <- hav2[apply(hav2 != 0, 1, all), ]
colnames(hav2)[1] <- "genes"

# Find columns ending with -AD and -C
ad_cols <- grep("-AD$", names(hav2), value = TRUE)
c_cols  <- grep("-C$", names(hav2), value = TRUE)

# Extract prefixes
ad_prefixes <- sub("-AD$", "", ad_cols)
c_prefixes  <- sub("-C$", "", c_cols)

# Find common prefixes
common_prefixes <- intersect(ad_prefixes, c_prefixes)

# new.data <- data.frame(row.names = hav2$genes)
# new.data[["Classical monocytes"]] <- hav2[["Classical.monocytes-AD"]] / hav2[["Classical.monocytes-C"]]

# Build result dataframe
result <- data.frame(row.names = hav2$genes)
for (prefix in common_prefixes) {
  ad_col <- paste0(prefix, "-AD")
  c_col  <- paste0(prefix, "-C")
  result[[prefix]] <- hav2[[ad_col]] / hav2[[c_col]]
}

result2 <- as.matrix(result)
result3 <- log2(result2)

data <- result3
# Define color breaks centered at 0, with more breaks near zero to enhance small differences
max_abs <- max(abs(data))
# Use more breaks between -1 and 1 for sensitivity
breaks <- c(seq(-max_abs, -1, length.out=30),
            seq(-1, 1, length.out=40),
            seq(1, max_abs, length.out=30))

# Create a blue-white-red color palette
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)
colors <- colorRampPalette(c("#003560", "white", "#7A1B00"))(length(breaks) - 1)

# Plot the heatmap
par(mar=c(16,5,2,2))
pheatmap(data,
         color = colors,
         breaks = breaks,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 7,
         # angle_col = "45",
         main = "AD-GWAS genes")


## Targets heatmap
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
myBreaks <- c(seq(min(result3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(result3)/paletteLength, max(result3), length.out=floor(paletteLength/2)))
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
pheatmap(result3, cluster_rows = T, show_rownames = T, 
                 show_colnames = T, #annotation_col = sampleTable2,
                 cluster_cols = T, annotation_names_col = T, annotation_legend = T,
                 color = colors,
                 main = "AD-GWAS genes", fontsize_row = 7,
                 #labels_row = as.expression(newnames),
                 # angle_col = 45,
                 breaks = myBreaks
                 # scale = "none"
) 

# hav3 <- hav2
# row.names(hav3) <- hav3$genes
# hav3 <- as.matrix(hav3[2:51])
# ## Targets heatmap
# paletteLength <- 500
# colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
# #png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
# heat <- pheatmap(hav3, cluster_rows = T, show_rownames = T, 
#                  show_colnames = T, #annotation_col = sampleTable2,
#                  cluster_cols = T, annotation_names_col = T, annotation_legend = T,
#                  color = colors,
#                  main = "AD-GWAS genes", fontsize_row = 7,
#                  #labels_row = as.expression(newnames),
#                  # angle_col = 45,
#                  scale = "row"
# ) 
# heat

############ Cytokine panel

# hav <- t(unique(read.delim("AD-GWAS.txt", header = F)))
# strips <- "^\\s+|\\s+$"
# hav <- gsub(strips,"",hav)
hav <- c("IL1A","IL1B","IL2","IL3","IL4","IL6","IL7","CXCL8","IL10","IL12A","IL15","IL16","IL17A","IL17F","IL18","IL22","IL23A","IL27",
         "IL32","IL33","IL34","IFNA1","IFNB1","IFNG","TGFB1","TNF")

Idents(prot.combined) <- "celltypes.disease"

avg.exp <- AverageExpression(
  prot.combined,
  assays = "RNA",
  features = hav,
  return.seurat = FALSE,
  group.by = "celltypes.disease",
  add.ident = NULL,
  layer = "data",
  # slot = deprecated(),
  verbose = TRUE
)

write.xlsx(as.data.frame(avg.exp$RNA), rowNames = T,file=paste0("cytokine_panel_avg.exp_PBMC.xlsx"))

hav2 <- read.xlsx("cytokine_panel_avg.exp_PBMC.xlsx")
hav2 <- hav2[apply(hav2 != 0, 1, all), ]
colnames(hav2)[1] <- "genes"

hav2b <- read.xlsx("cytokine_panel_avg.exp_PBMC.xlsx")
colnames(hav2b)[1] <- "genes"
row.names(hav2b) <- hav2b$genes
hav2 <- hav2b[2:51]
hav2 <- hav2+0.000000000001


# Find columns ending with -AD and -C
ad_cols <- grep("-AD$", names(hav2), value = TRUE)
c_cols  <- grep("-C$", names(hav2), value = TRUE)

# Extract prefixes
ad_prefixes <- sub("-AD$", "", ad_cols)
c_prefixes  <- sub("-C$", "", c_cols)

# Find common prefixes
common_prefixes <- intersect(ad_prefixes, c_prefixes)

# new.data <- data.frame(row.names = hav2$genes)
# new.data[["Classical monocytes"]] <- hav2[["Classical.monocytes-AD"]] / hav2[["Classical.monocytes-C"]]

## Build result dataframe
result <- data.frame(row.names = hav2$genes)
# result <- data.frame(row.names = row.names(hav2))
for (prefix in common_prefixes) {
  ad_col <- paste0(prefix, "-AD")
  c_col  <- paste0(prefix, "-C")
  result[[prefix]] <- hav2[[ad_col]] / hav2[[c_col]]
}

result2 <- as.matrix(result)
result3 <- log2(result2)

data <- result3
data[!is.finite(data)] <- 0
# Define color breaks centered at 0, with more breaks near zero to enhance small differences
max_abs <- max(abs(data))
# Use more breaks between -1 and 1 for sensitivity
breaks <- c(seq(-max_abs, -1, length.out=30),
            seq(-1, 1, length.out=40),
            seq(1, max_abs, length.out=30))

# Create a blue-white-red color palette
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)
colors <- colorRampPalette(c("#003560", "white", "#7A1B00"))(length(breaks) - 1)

# Plot the heatmap
par(mar=c(16,5,2,2))
pheatmap(data,
         color = colors,
         breaks = breaks,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 7,
         angle_col = "45",
         main = "AD-GWAS genes")


## Targets heatmap
result3[!is.finite(result3)] <- 0
result3 <- result3[,"Classical.monocytes"]
paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
myBreaks <- c(seq(min(result3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(result3)/paletteLength, max(result3), length.out=floor(paletteLength/2)))
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
par(mar=c(10,5,2,2))
pheatmap(result3, cluster_rows = F, show_rownames = T, 
         show_colnames = T, #annotation_col = sampleTable2,
         cluster_cols = F, annotation_names_col = T, annotation_legend = T,
         color = colors,
         main = "Cytokine panel", fontsize_row = 7,
         #labels_row = as.expression(newnames),
         angle_col = "45",
         breaks = myBreaks
         # scale = "none"
) 

# hav3 <- hav2
# row.names(hav3) <- hav3$genes
# hav3 <- as.matrix(hav3[2:51])
# ## Targets heatmap
# paletteLength <- 500
# colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
# #png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
# heat <- pheatmap(hav3, cluster_rows = T, show_rownames = T, 
#                  show_colnames = T, #annotation_col = sampleTable2,
#                  cluster_cols = T, annotation_names_col = T, annotation_legend = T,
#                  color = colors,
#                  main = "AD-GWAS genes", fontsize_row = 7,
#                  #labels_row = as.expression(newnames),
#                  # angle_col = 45,
#                  scale = "row"
# ) 
# heat

hav2 <- read.xlsx("cytokine_panel_avg.exp_PBMC.xlsx")
colnames(hav2)[1] <- "genes"
row.names(hav2) <- hav2$genes
result3 <- hav2[2:3]

paletteLength <- 500
colors <- colorRampPalette( rev(brewer.pal(11, "RdBu")) )(paletteLength)
myBreaks <- c(seq(min(result3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(result3)/paletteLength, max(result3), length.out=floor(paletteLength/2)))
#png(filename = "heatmap_DEGs_intercept.png", units = "px", res = 300, height = 3500, width = 3000)
pheatmap(result3, cluster_rows = T, show_rownames = T, 
         show_colnames = T, #annotation_col = sampleTable2,
         cluster_cols = F, annotation_names_col = T, annotation_legend = T,
         color = colors,
         main = "AD-GWAS genes", fontsize_row = 7,
         #labels_row = as.expression(newnames),
         # angle_col = 45,
         # breaks = myBreaks
         scale = "none"
) 

Idents(prot.combined) <- "celltypes"
VlnPlot(prot.combined, features = hav, idents = "Classical monocytes", group.by = "disease", flip = T, stack = T) +
  # scale_y_continuous(limits = c(0.000, 5)) + 
  theme(legend.position = "none")

DotPlot(prot.combined, features = rev(hav),
        # cols = "RdBu",
        cols = c("gray90","#9E1021"),
        col.max = 20,
        dot.scale = 10,
        # cluster.idents = T,
        scale = F,
        group.by = "disease",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()


###### Suppl to fig 2
## Violin plots

genes1 <- c("JAK1","IL15","IL15RA","JAK2","CD36","CREB1","STAT1","STAT2","STAT3","IRF1","IRF2","IRF3","TGFB1","TGFBR2","TGFBR1",
            "TNF","IFNGR1","IFNGR2","TLR2","TLR4","TLR6","MAPK14","PIAS1","PPARG","CD86","CD80",
            "SREBF1","SREBF2","RUNX3","PPARA","LRP1","CREBBP",
            "SPI1","CD14","JUN","JUNB","FOS","FOSL2","NFKB1","MYD88","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8") # Classical monocytes

genes1 <- c()

## genes compared by disease
my_comparisons <- list(c("AD","C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "disease",
                idents = "Classical monocytes"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.disease.pdf"),
         plot = p1, width = 6, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## genes compared by disease and sex

my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "sex.disease",
                idents = "Classical monocytes"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.sex.disease.pdf"),
         plot = p1, width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## genes compared by disease and sex and age

prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")

my_comparisons <- list(c("50-70_male_AD","50-70_male_C"),c("71-90_male_AD","71-90_male_C"),
                       c("50-70_female_AD","50-70_female_C"),c("71-90_female_AD","71-90_female_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "sex.age2.disease",
                idents = "Classical monocytes"
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

## genes compared by disease and age

my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "age.range.disease2",
                idents = "Classical monocytes"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_classic.mono_by.age.disease2.pdf"),
         plot = p1, width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

# #####  GO barplot
# 
# paths <- read.delim("reactome_classical.mono_upregulated_sel1.csv", header = T)
# colunas1 <- paths$Term[23:1]
# 
# plt <- ggplot(paths) + scale_x_continuous(
#   limits = c(0, 7),
#   breaks = seq(0, 7, by = 1), 
#   #expand = c(0, 0), # The horizontal axis does not extend to either side
#   position = "top"  # Labels are located on the top
# ) +
#   geom_col(aes(log10pval, factor(Term,levels=colunas1)), fill = paths$col, width = 0.8) + 
#   geom_text(
#     data = subset(paths, log10pval > 0),
#     aes(0, y = factor(Term,levels=colunas1), label = factor(Term,levels=colunas1)),
#     hjust = 0,
#     nudge_x = 0.1,
#     colour = "black",
#     #family = "Econ Sans Cnd",
#     size = 4
#   ) + geom_text(
#     data = subset(paths, log10pval < 0),
#     aes(0, y = Term, label = Term),
#     hjust = 1,
#     nudge_x = -0.1,
#     colour = "black",
#     #family = "Econ Sans Cnd",
#     size = 3
#   ) + xlab(expression("-Log"[10]*"("*italic("p")*"-value)")) +
#   ggtitle("Effector CD8+ T cells upregulated genes") + ylab("Term") + #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
#   theme(
#     ## Set background color to white
#     panel.background = element_rect(fill = "white"),
#     ## Set the color and the width of the grid lines for the horizontal axis
#     panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
#     ## Remove tick marks by setting their length to 0
#     axis.ticks.length = unit(0, "mm"),
#     ## Remove labels from the vertical axis
#     axis.text.y = element_blank(),
#   )
# plt
# ggsave(filename = "Effector.CD8.T_up_DEGs_GOBP_plot.pdf",
#        plot = plt,
#        width = 12,
#        height = 4,
#        device = "pdf")

## Featureplots
FeaturePlot(prot.combined, raster = F, pt.size = 0.05, features = "IFNG", 
            split.by = "disease", 
            reduction = "umap")

### Dotplots
genes1 <- c("IFNGR2","JAK1","JAK2","STAT1","STAT2","IRF2",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K4","MAP4K5","MAPK1","MAPK8",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # nonclassic Monocytes

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Nonclassical monocytes",
) + RotatedAxis() + coord_flip()

genes1 <- c(#"IFNGR1",
  "JAK1","JAK2","IRF2",
  "PIAS1","RUNX3","PPARG","CD36","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
  "MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K4","MAP3K5", 
  "MAP3K7","MAP3K20","MAP4K3","MAP4K5","MAPK8",
  "CREB1",
  "IL15"#,"CD86"
) # intermediate Monocytes

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "age.range.disease2",
        idents = "Intermediate monocytes",
        # Classical monocytes
) + RotatedAxis() + coord_flip()


prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")
my_comparisons <- c("50-70_male_C","50-70_male_AD",
                    "50-70_female_C","50-70_female_AD",
                    "71-90_male_C","71-90_male_AD",
                    "71-90_female_C","71-90_female_AD")
prot.combined$sex.age2.disease <- factor(prot.combined$sex.age2.disease, levels = my_comparisons)

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "sex.age2.disease",
        idents = "Intermediate monocytes",
) + RotatedAxis() + coord_flip()


############# Suppl to fig 3
######## Fig 3B alternative

prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")
my_comparisons <- c("50-70_male_C","50-70_male_AD",
                    "50-70_female_C","50-70_female_AD",
                    "71-90_male_C","71-90_male_AD",
                    "71-90_female_C","71-90_female_AD")
prot.combined$sex.age2.disease <- factor(prot.combined$sex.age2.disease, levels = my_comparisons)

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "sex.age2.disease",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()

######## Violin plots

genes1 <- c("PRF1","GZMA","GZMH","GZMK","GZMB","TOX","NKG7","IL7R","IRF2","IRF1","IRF3","IL2RG","IL2RB","STAT4","STAT1","STAT2","STAT3","JAK1","JAK2",
            "MAPK1","CD8A","CD8B","TRAC","TRBC1","CD28") # Effector CD8

## genes compared by disease
my_comparisons <- list(c("AD","C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "disease",
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+0.5))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.disease.pdf"),
         plot = p1, width = 6, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## genes compared by disease and sex

my_comparisons <- list(c("male_AD","male_C"),c("female_AD","female_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "sex.disease",
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.sex.disease.pdf"),
         plot = p1, width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## genes compared by disease and sex and age

prot.combined$sex.age2.disease <- paste(prot.combined$age.range2,prot.combined$sex,prot.combined$disease,sep = "_")

my_comparisons <- list(c("50-70_male_AD","50-70_male_C"),c("71-90_male_AD","71-90_male_C"),
                       c("50-70_female_AD","50-70_female_C"),c("71-90_female_AD","71-90_female_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "sex.age2.disease",
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+2))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.sex.age2.disease.pdf"),
         plot = p1, width = 12, height = 10,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

## genes compared by disease and age

my_comparisons <- list(c("50-70_AD","50-70_C"),c("71-90_AD","71-90_C"))

for (i in 1:length(genes1)){
  p1 <- VlnPlot(prot.combined, features = genes1[i],
                pt.size = 0.05, raster = F, group.by = "age.range.disease2",
                idents = "Effector CD8+ T cells"
  ) 
  p1 <- p1 + stat_summary(fun = mean, geom='point', size = 35, colour = "black", shape = 95) +
    scale_y_continuous(limits = c(0.000, max(p1[[1]][["data"]][[genes1[i]]])+1))+
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
  p1$layers[[2]]$aes_params$alpha <- 0.3
  ggsave(filename = paste0("vln_stats_",genes1[i],"_mean_cd8.effect_by.age.disease2.pdf"),
         plot = p1, width = 8, height = 8,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}


#####  GO barplot

paths <- read.delim("GO_BP_effector_t.cd8_upregulated_sel.csv", header = T)
colunas1 <- paths$Term[23:1]

plt <- ggplot(paths) + scale_x_continuous(
  limits = c(0, 7),
  breaks = seq(0, 7, by = 1), 
  #expand = c(0, 0), # The horizontal axis does not extend to either side
  position = "top"  # Labels are located on the top
) +
  geom_col(aes(log10pval, factor(Term,levels=colunas1)), fill = paths$col, width = 0.8) + 
  geom_text(
    data = subset(paths, log10pval > 0),
    aes(0, y = factor(Term,levels=colunas1), label = factor(Term,levels=colunas1)),
    hjust = 0,
    nudge_x = 0.1,
    colour = "black",
    #family = "Econ Sans Cnd",
    size = 4
  ) + geom_text(
    data = subset(paths, log10pval < 0),
    aes(0, y = Term, label = Term),
    hjust = 1,
    nudge_x = -0.1,
    colour = "black",
    #family = "Econ Sans Cnd",
    size = 3
  ) + xlab(expression("-Log"[10]*"("*italic("p")*"-value)")) +
  ggtitle("Effector CD8+ T cells upregulated genes") + ylab("Term") + #scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    ## Set background color to white
    panel.background = element_rect(fill = "white"),
    ## Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = "#A8BAC4", linewidth = 0.3),
    ## Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
  )
plt
ggsave(filename = "Effector.CD8.T_up_DEGs_GOBP_plot.pdf",
       plot = plt,
       width = 12,
       height = 4,
       device = "pdf")


####### Other supplementaries

########### signature by cohort

Idents(prot.combined) <- "celltypes"
my_comparisons <- list(c("Harvard_AD","Harvard_C"),c("Sichuan_AD","Sichuan_C"),
                       c("Chicago_AD","Chicago_C"))
colors <- brewer.pal(n=3,name = "Dark2")
colors <- colors[c(2,1,2,1,2,1)]
order1 <- c("Harvard_C","Harvard_AD","Sichuan_C","Sichuan_AD","Chicago_C","Chicago_AD")
prot.combined$cohort.disease <- factor(prot.combined$cohort.disease, levels = order1)

gene1 <- "IL15"
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors, 
              idents = c("Classical monocytes")
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = mean, geom='point', size = 25, colour = "black", shape = 95) +
  scale_y_continuous(limits = c(0.001, max(p1[[1]][["data"]][[gene1]])+0.5))+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1

genes1 <- c("IFNGR1","IFNGR2","JAK1","JAK2","STAT1","STAT2","STAT3","IRF2",#"IRF3","IRF1",
            "PIAS1","RUNX3","PPARA","PPARG","CD36","TLR2","TLR4","MYD88",#"TLR6","LRP1","SREBF1","SREBF2",
            "MAPK14","MAP2K1","MAP2K4","MAP2K5","MAP2K6","MAP3K1","MAP3K3","MAP3K4","MAP3K5", 
            "MAP3K7","MAP3K14","MAP3K20","MAP4K3","MAP4K5","MAPK1","MAPK8", #"MAP4K4",
            "CREB1","CREBBP",
            "IL15","IL15RA"#,"CD86"
) # classic Monocytes

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "cohort.disease",
        idents = "Classical monocytes",
) + RotatedAxis() + coord_flip()

prot.combined <- AddModuleScore(
  prot.combined,
  features = genes1,
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
# colors <- brewer.pal(n=6,name = "Dark2")
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors,
              idents = "Classical monocytes"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.0000, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1

genes1 <- c( "CD8B","CD8A", # "IL2RA","TNF","STAT4",
  "TIGIT","LTA","SLAMF6","IL2RB","IL2RG","EOMES", #"IRF2","GNLY", #"IRF3","IRF4", #"LTB", "IRF1",
  "PRF1","GZMH","GZMA", "GZMM", # "GZMK","GZMB",
  "NKG7","TYROBP","KLRB1",
  "BATF","TOX","IFNG","CD244","TBX21" #,
  #"ENTPD1","FASLG",
  #"ZAP70","LCK","TCF7","CD28","NFATC2" #,"NR4A1","NFATC1","FYN",
  
) # effector cd8

DotPlot(prot.combined, features = rev(genes1),
        cols = "RdBu",
        col.max = 20,
        dot.scale = 15,
        # cluster.idents = T,
        #scale = F,
        group.by = "cohort.disease",
        idents = "Effector CD8+ T cells",
) + RotatedAxis() + coord_flip()

prot.combined <- AddModuleScore(
  prot.combined,
  features = genes1,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = "RNA",
  name = "Cytotoxicity.signature",
  seed = 1,
  search = FALSE,
  slot = "data"
)

gene1 <- "Cytotoxicity.signature1"
# colors <- brewer.pal(n=6,name = "Dark2")
p1 <- VlnPlot(prot.combined, features = gene1,
              pt.size = 0.05, raster = F, group.by = "cohort.disease", cols = colors,
              idents = "Effector CD8+ T cells"
) + theme(legend.position = "none")
p1 <- p1 + stat_summary(fun = median, geom='point', size = 15, colour = "black", shape = "\U2014") +
  #stat_summary(fun = mean, geom='point', size = 5, colour = "darkred") + #geom_boxplot(width=0.1) +
  scale_y_continuous(limits = c(0.001, max(p1[[1]][["data"]][[gene1]])+0.0*max(p1[[1]][["data"]][[gene1]]))) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") #+ # Add pairwise comparisons p-value
p1$layers[[2]]$aes_params$alpha <- 0.01
p1


########## clusters frequencies

df_freqs <- read.csv("freq.num.clusters_disease_patient_prot.combined.csv", header = T)
celltypes2 <- colnames(df_freqs)[11:93]
coluna1 <- c("AD","C")
# coluna1 <- c("AD_male","C_male","AD_female","C_female")
# coluna1 <- c("AD_50-70","C_50-70","AD_71-90","C_71-90")
# coluna1 <- c("AD_50-70_male","C_50-70_male","AD_50-70_female","C_50-70_female",
#              "AD_71-90_male","C_71-90_male","AD_71-90_female","C_71-90_female")

# my_comparisons <- list(c("AD_50-70_male","C_50-70_male"),c("AD_50-70_female","C_50-70_female"),
#                        c("AD_71-90_male","C_71-90_male"),c("AD_71-90_female","C_71-90_female"))
#my_comparisons <- list(c("AD_male","C_male"),c("AD_female","C_female"),c("AD_female","AD_male"),c("C_female","C_male"))
#my_comparisons <- list(c("AD_71-90","C_71-90"),c("AD_50-70","C_50-70"),c("AD_50-70","AD_71-90"),c("C_71-90","C_50-70"))
my_comparisons <- list(c("AD","C"))


for (i in 1:length(celltypes2)){
  p1 <- ggplot(df_freqs, aes_string(x="factor(disease,levels=coluna1)", y=celltypes2[i], color="disease")) + 
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
  ggsave(filename = paste0("vln_stats_",celltypes2[i],"_by.disease.pdf"),
         plot = p1,
         width = 4,
         height = 6,
         #dpi = 600,
         device = "pdf")
  rm(p1)
}

########## cell types frequencies

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


######## DEGs by clusters

clusters <- as.character(t(read.delim("cluster.order2.txt", header = F)))
clusters_cells <- as.character(t(read.delim("clusters.order2_celltypes2.txt", header = F)))

data <- data.frame(matrix(ncol = 5,nrow = 0))

for (i in 1:length(clusters)){
  list1 <- read.xlsx(paste0("DEGs/age2.sex.cohort.regress/clusters/MAST_clusters_lat.vars.age2.sex.cohort.Harvard_pct0.3_ADxHC_clus",clusters[i],"_DEGs.xlsx"))
  data[i,1] <- clusters[i]
  data[i,2] <- length(row.names(list1))
  data[i,3] <- length(list1$avg_log2FC[list1$avg_log2FC > 0])
  data[i,4] <- length(list1$avg_log2FC[list1$avg_log2FC < 0])
  data[i,5] <- clusters_cells[i]
  rm(list1)
}

colnames(data) <- c("clusters_num","DEGs","upregulated","downregulated","clusters")

coluna1 <- data$clusters
coluna1 <- rev(coluna1) 
data$clusters <- factor(data$clusters, levels=coluna1)

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
    axis.text.y = element_blank(),
  )

p2 <- ggplot(data, aes(y=factor(clusters,levels=coluna1), x=downregulated)) +
  geom_col(fill = "skyblue") + scale_x_continuous(
    #limits = c(0, 5000),
    #breaks = seq(0, 5000, by = 500),
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
    ## Remove labels from the vertical axis
    axis.text.y = element_blank(),
  )

ggarrange(p1, p2)


















