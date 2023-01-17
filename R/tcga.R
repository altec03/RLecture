library(tidyverse)

load(here::here("RData/tnbcType.RData"))

library(maftools)

tcgaMaf <- read.maf(here::here("data/data_mutations.txt"))

mafTumorSampleBarcode <- gsub("-[0-9]{2}", "", tcgaMaf@clinical.data$Tumor_Sample_Barcode)

tnbcType <- tnbcType[tnbcType$Tumor_Sample_Barcode %in% mafTumorSampleBarcode,]

tcgaMaf@data$Tumor_Sample_Barcode <- gsub("-[A-Z0-9]{2}$", '', tcgaMaf@data$Tumor_Sample_Barcode)

tnbcMaf <- subsetMaf(maf = tcgaMaf, tsb = tnbcType$Tumor_Sample_Barcode[!is.na(tnbcType$TIL_Group)], mafObj = TRUE, isTCGA = FALSE)

tnbcType$Tumor_Sample_Barcode[!is.na(tnbcType$TIL_Group)]

genes <- "TP53 MCL1 MYC CCNE1 CDKN2A RB1 BRCA1 BRCA2 ATM FANCD2 PIK3CA PTEN ERBB2 EGFR NF1 NOTCH1 NOTCH2 NOTCH3 FBXW7 BRAF HRAS MAP2K1 KRAS MAPK9 DNMT1 ASH2L ASXL2 ASXL3 EZH2 ARID1A ARID1B KDM6A BAP1 SMARCA4 SMARCAD1 BCL7C HCFC1 B2M CD274 PDCD1LG2 CTLA4 IDO1"
genes <- str_split(genes, " ")|> unlist()

cnv <- read_delim(here::here("data/data_cna.txt"))

cnTable <- cnv |>
  pivot_longer(cols = -1:-2, names_to = "Sample_name", values_to = "CN")|>
  filter(Hugo_Symbol %in% genes)|>
  mutate(Sample_name= gsub("-[0-9]{2}", "", Sample_name))|>
  filter(CN != 0,
         Sample_name %in% tnbcType$Tumor_Sample_Barcode)|>
  mutate(CN = case_when(CN ==-2 ~ 'DeepDel',
                        CN == -1 ~ 'Del',
                        CN == 1 ~ 'ShallowAmp',
                        CN == 2 ~ 'Amp'
  ))|>
  select(Hugo_Symbol, Sample_name, CN)|>
  filter(Sample_name %in% tnbcType$Tumor_Sample_Barcode[!is.na(tnbcType$TIL_Group)])

tnbcMaf <- read.maf(maf = tnbcMaf@data,
                    cnTable = cnTable,
                    clinicalData = tnbcType)

save(tnbcMaf, file = here::here("RData/tnbcMaf.RData"))

