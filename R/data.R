mutation <- readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))
clinical <- readr::read_csv(here::here("data/Pathology_NGS_1st.csv"))

library(tidyverse)
colnames(clinical)

clinical <- clinical |>
  select(-name, -patientID, -rel_pathology_num, -starts_with(c("tsv", "img")),
         -sendEMR, -sendEMRDate, -examin, -recheck, -appoint_doc, -key_block) 
