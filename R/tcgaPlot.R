library(plotgardener)
library(ggpubr)
library(tidyverse)

load(here::here("RData/tnbcMafFirehose.RData"))
load(here::here("RData/tnbcMaf.RData"))
load(here::here("RData/pal.RData"))

genes <- "TP53 MCL1 MYC CCNE1 CDKN2A RB1 BRCA1 BRCA2 ATM FANCD2 PIK3CA PTEN ERBB2 EGFR NF1 NOTCH1 NOTCH2 NOTCH3 FBXW7 BRAF HRAS MAP2K1 KRAS MAPK9 DNMT1 ASH2L ASXL2 ASXL3 EZH2 ARID1A ARID1B KDM6A BAP1 SMARCA4 SMARCAD1 BCL7C HCFC1 B2M CD274 PDCD1LG2 CTLA4 IDO1"
genes <- str_split(genes, " ")|> unlist()
pathways = data.frame(
  genes,
  Pathway = rep(c(
    "Cell cycle", "DNA repair", "PI3K", "NOTCH", "MAPK", "Epigenetic", "Immune"
    ),
    c(6, 4, 5, 4, 5, 13, 5)),
    stringsAsFactors = FALSE
)

anno_color <- myPal
names(anno_color) <- c("LAR", "LP", "LI", "LD")
list_anno_color <- list(TIL_Group = anno_color)
oncoplot(tnbcMaf, 
         clinicalFeatures = c("TIL_Group"),
         sortByAnnotation = TRUE,
         groupAnnotationBySize = FALSE,
         annotationOrder = c("LAR", "LP", "LI", "LD"),
         topBarData = "HRD",
         drawColBar = FALSE,
         includeColBarCN = FALSE,
         bgCol = "white",
         borderCol = "white",
         gene_mar = 6,
         legend_height = 4,
         anno_height = 0.2,
         fontSize = 0.6,
         titleFontSize = 1.0,
         legendFontSize = 0.8,
         annotationFontSize = 0.8,
         annotationColor = list_anno_color,
         pathways = pathways)

my_comparisons <- list( c("LAR", "LP"), c("LAR", "LI"), c("LAR", "LD") )

tnbcType|>
  filter(!is.na(TIL_Group))|>
  ggpubr::ggviolin(x = "TIL_Group", y = "HRD", 
                   fill = "TIL_Group",
                   palette = myPal,
                   add = "boxplot",
                   alpha = 0.5,
                   xlab = "", ylab = "Genomic scar signature scores",
                   ylim = c(0,170),
                   font.label = list(size = 4),
                   remove = "legend"
                  )+
  ggpubr::stat_compare_means(comparisons = my_comparisons,
                             label.y = c(100, 115, 130),
                             label = "p.signif",
                             method = "wilcox.test",
                             vjust = 0.2)+
  ggpubr::stat_compare_means(label.y = 155, label.x = 1.5) +
  rremove("legend")

