library(tidyverse)
load(here::here("RData/tnbcCnvFirehose.RData"))
load(here::here("RData/tnbcType.RData"))
load(here::here("RData/pal.RData"))

genes <- "TP53 MCL1 MYC CCNE1 CDKN2A RB1 BRCA1 BRCA2 ATM FANCD2 PIK3CA PTEN ERBB2 EGFR NF1 NOTCH1 NOTCH2 NOTCH3 FBXW7 BRAF HRAS MAP2K1 KRAS MAPK9 DNMT1 ASH2L ASXL2 ASXL3 EZH2 ARID1A ARID1B KDM6A BAP1 SMARCA4 SMARCAD1 BCL7C HCFC1 B2M CD274 PDCD1LG2 CTLA4 IDO1"
genes <- str_split(genes, " ")|> unlist()

dataChiCna <- cnv |>
  pivot_longer(cols = -1:-2, names_to = "Tumor_Sample_Barcode", values_to = "CN")|>
  mutate(Tumor_Sample_Barcode= gsub("-[0-9]{2}", "", Tumor_Sample_Barcode))|>
  filter(Tumor_Sample_Barcode %in% tnbcType$Tumor_Sample_Barcode)|>
  filter(Hugo_Symbol %in% genes)

dataDel <- dataChiCna|>
  mutate(CN = case_when(CN ==-2 ~ 'Del',
                        CN == -1 ~ 'Del',
                        CN == 0 ~ "Not Del",
                        CN == 1 ~ "Not Del",
                        CN == 2 ~ "Not Del"
  ))|>
  select(Hugo_Symbol, Tumor_Sample_Barcode, CN)|>
  filter(Tumor_Sample_Barcode %in% tnbcType$Tumor_Sample_Barcode[!is.na(tnbcType$TIL_Group)])


dataFisherDel <- dataDel |>
  left_join(tnbcType |>
              select(Tumor_Sample_Barcode, TIL_Group))

fisherDel <- dataFisherDel |>
  group_by(Hugo_Symbol) |>
  nest()|>
  mutate(delNA = map(data, ~ sum(.x$CN=="Del")))|>
  filter(delNA > 20)|>
  mutate(fisher = map(data, ~ fisher.test(.x$CN, .x$TIL_Group)))|>
  mutate(tidyfisher = map(fisher, broom::tidy)) |>
  unnest(tidyfisher) 


####################### AMP ########################################

dataAmp <- dataChiCna|>
  mutate(CN = case_when(CN ==-2 ~ 'Not Amp',
                        CN == -1 ~ 'Not Amp',
                        CN == 0 ~ "Not Amp",
                        CN == 1 ~ "Not Amp",
                        CN == 2 ~ "Amp"
  ))|>
  select(Hugo_Symbol, Tumor_Sample_Barcode, CN)|>
  filter(Tumor_Sample_Barcode %in% tnbcType$Tumor_Sample_Barcode[!is.na(tnbcType$TIL_Group)])


dataFisherAmp <- dataAmp |>
  left_join(tnbcType |>
              select(Tumor_Sample_Barcode, TIL_Group))

tempData <- dataFisherAmp |>
  group_by(Hugo_Symbol) |>
  nest()

fisherAmp <- dataFisherAmp |>
  group_by(Hugo_Symbol) |>
  nest()|>
  mutate(ampNA = map(data, ~ sum(.x$CN=="Amp")))|>
  filter(ampNA > 20)|>
  mutate(fisher = map(data, ~ fisher.test(.x$CN, .x$TIL_Group)))|>
  mutate(tidyfisher = map(fisher, broom::tidy)) |>
  unnest(tidyfisher) 

##################### Mutation #########################3

library(tidyverse)
load(here::here("RData/tnbcMafFirehose.RData"))
mafData <- tnbcMaf@data[,-79:-80]
mafData <- mafData|>
  filter(Variant_Type != "CNV")

PIK3CA <- tnbcType |>
  left_join(mafData) |>
  select(Tumor_Sample_Barcode, TIL_Group, Hugo_Symbol)|>
  group_by(Tumor_Sample_Barcode)|>
  summarise(PIK3CA =sum(Hugo_Symbol == "PIK3CA"))|>
  mutate(PIK3CA = PIK3CA ==0,
         PIK3CA = factor(PIK3CA,
                         levels = c(TRUE, FALSE),
                         labels = c("Wild-type", "Mutant")))|>
  left_join(tnbcType)|>
  select(PIK3CA, TIL_Group)|>
  table()

fisherPik3ca <- PIK3CA|>
  fisher.test()|>
  broom::tidy()|>
  mutate(Hugo_Symbol = "PIK3CA")
  
PTEN <- tnbcType |>
  left_join(mafData) |>
  select(Tumor_Sample_Barcode, TIL_Group, Hugo_Symbol)|>
  group_by(Tumor_Sample_Barcode)|>
  summarise(PTEN =sum(Hugo_Symbol == "PTEN"))|>
  mutate(PTEN = PTEN ==0,
         PTEN = factor(PTEN,
                         levels = c(TRUE, FALSE),
                         labels = c("Wild-type", "Mutant")))|>
  left_join(tnbcType)|>
  select(PTEN, TIL_Group)|>
  table()


fisherPTEN <- PTEN|>
  fisher.test()|>
  broom::tidy()|>
  mutate(Hugo_Symbol = "PTEN")

tp53 <- tnbcType |>
  left_join(mafData) |>
  select(Tumor_Sample_Barcode, TIL_Group, Hugo_Symbol)|>
  group_by(Tumor_Sample_Barcode)|>
  summarise(TP53 =sum(Hugo_Symbol == "TP53"))|>
  mutate(TP53 = TP53 ==0,
         TP53 = factor(TP53,
                       levels = c(TRUE, FALSE),
                       labels = c("Wild-type", "Mutant")))|>
  left_join(tnbcType)|>
  select(TP53, TIL_Group)|>
  table()


fisherTp53 <- tp53|>
  fisher.test()|>
  broom::tidy()|>
  mutate(Hugo_Symbol = "TP53")
  

################### FDR #########################################

genomeComparison <- fisherAmp |>
  bind_rows(fisherDel) |>
  bind_rows(fisherPik3ca)|>
  bind_rows(fisherPTEN)|>
  bind_rows(fisherTp53)|>
  mutate(fdr = p.adjust(p.value, method = "fdr"))

genomeComparison$fdr <- p.adjust(genomeComparison$p.value, method = "fdr")

genomeComparison|>
  filter(fdr <= 0.05)|>
  select(Hugo_Symbol, p.value, fdr)


######################## bar plots #####################################


pik3caPlot <- PIK3CA|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(PIK3CA == "Mutant")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("PIK3CA mutation")))+
  xlab("")

mycPlot <- dataFisherAmp |>
  filter(Hugo_Symbol == "MYC") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Amp")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("MYC amplification")))+
  xlab("")


notch2Plot <- dataFisherAmp |>
  filter(Hugo_Symbol == "NOTCH2") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Amp")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("NOTCH2 \namplification")))+
  xlab("")

PDCD1LG2Plot <- dataFisherDel |>
  filter(Hugo_Symbol == "PDCD1LG2") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Del")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("PDCD1LG2(PD-L2)\n deletion")))+
  xlab("")

bap1Plot <- dataFisherDel |>
  filter(Hugo_Symbol == "BAP1") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Del")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("BAP1 deletion")))+
  xlab("")

cd274Plot <- dataFisherDel |>
  filter(Hugo_Symbol == "CD274") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Del")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("CD274(PD-L1)\n deletion")))+
  xlab("")

B2MPlot <- dataFisherDel |>
  filter(Hugo_Symbol == "B2M") |>
  select(CN, TIL_Group) |>
  table()|>
  prop.table(margin = 2)|>
  round(3)|>
  as_tibble()|>
  filter(CN == "Del")|>
  ggpubr::ggbarplot(x="TIL_Group", y="n",
                    fill = "TIL_Group",
                    label = FALSE, 
                    palette = myPal,
                    position = position_dodge(0.9))+
  theme(legend.position = "none")+
  ylab(expression(paste("B2M deletion")))+
  xlab("")

save(B2MPlot, bap1Plot, cd274Plot, mycPlot, notch2Plot, PDCD1LG2Plot, pik3caPlot,
     file = here::here("RData/barplotFirehose.RData"))
