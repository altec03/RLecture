

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SummarizedExperiment")
BiocManager::install("sesame")

library(TCGAbiolinks)
library(SummarizedExperiment)

#-----------------------------------
# STEP 1: Search, download, prepare |
#-----------------------------------
# 1.1 - DNA methylation
# ----------------------------------
query.met <- GDCquery(
  project = "TCGA-ACC", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

GDCdownload(
  query = query.met, 
  files.per.chunk = 20,
  directory = "case3/GDCdata"
)

acc.met <- GDCprepare(
  query = query.met,
  save = FALSE,
  directory = "case3/GDCdata"
)

#-----------------------------------
# 1.2 - RNA expression
# ----------------------------------
query.exp <- GDCquery(
  project = "TCGA-ACC", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

GDCdownload( 
  query = query.exp,
  files.per.chunk = 20,
  directory = "case3/GDCdata"
)

acc.exp <- GDCprepare(
  query = query.exp, 
  save = FALSE,
  directory = "case3/GDCdata"
)

# na.omit
acc.met <- acc.met[rowSums(is.na(assay(acc.met))) == 0,]

# Volcano plot
acc.met <- TCGAanalyze_DMC(
  data = acc.met, 
  groupCol = "paper_MethyLevel",#
  group1 = "CIMP-high",
  group2="CIMP-low",
  p.cut = 10^-5,
  diffmean.cut = 0.25,
  legend = "State",
  plot.filename = "case3/CIMP-highvsCIMP-low_metvolcano.png"
)

acc.exp.aux <- subset(
  acc.exp, 
  select = colData(acc.exp)$paper_MethyLevel %in% c("CIMP-high","CIMP-low")
)

idx <- colData(acc.exp.aux)$paper_MethyLevel %in% c("CIMP-high")
idx2 <- colData(acc.exp.aux)$paper_MethyLevel %in% c("CIMP-low")

dataPrep <- TCGAanalyze_Preprocessing(
  object = acc.exp.aux, 
  cor.cut = 0.6
)

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  qnt.cut = 0.25,
  method = 'quantile'
)

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,idx],
  mat2 = dataFilt[,idx2],
  Cond1type = "CIMP-high",
  Cond2type = "CIMP-low",
  method = "glmLRT"
)

TCGAVisualize_volcano(
  x = dataDEGs$logFC,
  y = dataDEGs$FDR,
  filename = "case3/Case3_volcanoexp.png",
  x.cut = 3,
  y.cut = 10^-5,
  names = rownames(dataDEGs),
  color = c("black","red","darkgreen"),
  names.size = 2,
  xlab = " Gene expression fold change (Log2)",
  legend = "State",
  title = "Volcano plot (CIMP-high vs CIMP-low)",
  width = 10
)

#------------------------------------------
# 2.4 - Starburst plot
# -----------------------------------------
# If true the argument names of the genes in circle 
# (biologically significant genes, has a change in gene
# expression and DNA methylation and respects all the thresholds)
# will be shown
# these genes are returned by the function see starburst object after the function is executed
starburst <- TCGAvisualize_starburst(
  met = acc.met, 
  exp = dataDEGs,
  genome = "hg19",
  group1 = "CIMP-high",
  group2 = "CIMP-low",
  filename = "case3/starburst.png",
  met.platform = "Illumina Human Methylation 450",
  met.p.cut = 10^-5,
  exp.p.cut = 10^-5,
  diffmean.cut = 0.25,
  logFC.cut = 3,
  names = FALSE, 
  height = 10,
  width = 15,
  dpi = 300
)

