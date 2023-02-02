library(tidyverse)

clinical <-
  readr::read_csv(here::here("data/Pathology_NGS_1st.csv"))

colnames(clinical)

class(clinical$prescription_date)
class(clinical$createDate)

clinical |>
  mutate(prescription_date = lubridate::ymd(prescription_date)) |>
  select(prescription_date, id, name) |>
  view()

library(lubridate)


clinical <- clinical |>
  #select(prescription_date, createDate)|>
  mutate(prescription_date = lubridate::ymd(prescription_date)) |>
  mutate(
    presc_year = year(prescription_date),
    presc_month = month(prescription_date),
    presc_week = week(prescription_date),
    year_month = floor_date(prescription_date, "month")
  )


clinical |>
  ggplot2::ggplot(aes(x = year_month)) +
  geom_bar()

clinical |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x = year_month, color = "blue")) +
  geom_bar(color = "blue") +
  geom_line(stat = "count")#, color = "red")+
xlab("Year (Monthly)") +
  ylab("Number of NGS test \n (ThermoFisher)") +
  theme_minimal()

clinical |>
  filter(year_month > "2020-10-01") |>
  group_by(year_month) |>
  summarise(n()) |>
  ggplot(aes(x = year_month, y = `n()`)) +
  geom_line()

clinical |> view()

class(clinical$pathological_dx)
summary(clinical$pathological_dx)

clinicalTumorType <- clinical |>
  mutate(
    diagnosis = as.factor(pathological_dx),
    cancer_type = str_to_lower(tumor_type),
    cancer_type = as.factor(cancer_type),
    cancer_type = case_when(
      str_detect(cancer_type, "colo|rectal") ~ "Colorectal Cancer",
      str_detect(cancer_type, "lung") ~ "Lung Cancer",
      str_detect(cancer_type, "glioblastoma|glioma|ependymoma|astrocytoma") ~ "Brain Cancer",
      str_detect(cancer_type, "pancrea") ~ "Pancreas Cancer",
      str_detect(cancer_type, "gastric") ~ "Stomach Cancer",
      str_detect(cancer_type, "breast") ~ "Breast Cancer",
      is.na(cancer_type) ~ NA_character_,
      TRUE ~ "Other Cancer"
    ),
    cancer_type = as.factor(cancer_type),
    cancer_type = fct_infreq(cancer_type)
  )

clinicalTumorType |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x = year_month)) + #, color = cancer_type) +
  #geom_bar()+
  geom_line(stat = "count") +
  facet_wrap(vars(cancer_type), ncol = 3) + #, scales = "free_y")+
  xlab("Year (Monthly)") +
  ylab("Number of NGS test \n (ThermoFisher)") +
  theme_minimal()

mutation <-
  readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))

data <- mutation |>
  full_join(clinicalTumorType, by = "pathology_num")

#######################################################
library(tidyverse)
load(here::here("RData/data.RData"))
skimr::skim(data)

data |>
  mutate(gene = as.factor(gene),
         gene = fct_lump(gene, 10),
         gene = fct_infreq(gene)) |>
  filter(gene != "Other") |>
  ggplot(aes(x = gene)) +
  geom_bar() +
  coord_flip() +
  facet_wrap( ~ cancer_type)

data |>
  ggplot(aes(msiscore)) +
  geom_histogram()

data |>
  mutate(
    msiscore_number = parse_number(msiscore),
    msi_status = str_extract(msiscore, "MS[SI]"),
    msi_status = factor(msi_status)
  ) |>
  select(msiscore, msi_status, msiscore_number) |>
  skimr::skim()

data <- data |>
  mutate(
    msiscore_number = parse_number(msiscore),
    msi_status = str_extract(msiscore, "MS[SI]"),
    msi_status = factor(msi_status)
  )


data |>
  ggplot(aes(msi_status, msiscore_number)) +
  geom_violin()

data |>
  filter(!(msi_status == "MSS" & msiscore_number > 200)) |>
  group_by(msi_status) |>
  summarise(
    n = n(),
    min_msiScore = min(msiscore_number),
    max_msiScore = max(msiscore_number),
    min_TMB = min(tumorburden, na.rm = TRUE),
    max_TMB = max(tumorburden, na.rm = TRUE)
  )

data |>
  filter(cancer_type == "Colorectal Cancer",
         is.na(gene)) |> view()

data |>
  group_by(gene) |>
  nest()

dfNested <- data |>
  group_by(pathology_num) |>
  distinct(gene, .keep_all = TRUE) |>
  ungroup() |>
  pivot_wider(
    id_cols = pathology_num,
    names_from = gene,
    values_from = nucleotide_change,
    values_fill = "wild_type"
  ) |>
  select(-"NA")|>
  left_join(dataTidy, by = "pathology_num")|>
  #colnames()
  pivot_longer(cols = USP9X:PMS1,
               names_to = "gene",
               values_to = "nucleotide_change") |>
  group_by(gene)|>
  nest()

dataTidy <- data |>
  select(pathology_num, cancer_type, starts_with("msi"), tumorburden)|>
  distinct_all()

dataNested <- dataTidy |>
  full_join(dfNested, by = "pathology_num")

save(dfNested, file = here::here("RData/dfNested.RData"))

data |>
  filter(!is.na(pathology_num),!is.na(gene)) |>
  group_by(pathology_num) |>
  distinct(gene, .keep_all = TRUE) |>
  ungroup() |>
  pivot_wider(
    id_cols = gene,
    names_from = pathology_num,
    values_from = nucleotide_change,
    values_fill = "wild_type"
  )

view()