library(tidyverse)

clinical <- readr::read_csv(here::here("data/Pathology_NGS_1st.csv"))

colnames(clinical)

class(clinical$prescription_date)
class(clinical$createDate)

clinical |>
  mutate(prescription_date = lubridate::ymd(prescription_date))|>
  select(prescription_date, id, name)|>
  view()

library(lubridate)


clinical <- clinical |>
  #select(prescription_date, createDate)|>
  mutate(prescription_date = lubridate::ymd(prescription_date)) |>
  mutate(presc_year = year(prescription_date),
         presc_month = month(prescription_date),
         presc_week = week(prescription_date),
         year_month = floor_date(prescription_date, "month"))


clinical |>
  ggplot2::ggplot(aes(x=year_month)) +
  geom_bar()

clinical |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x=year_month, color = "blue")) +
  geom_bar(color = "blue")+
  geom_line(stat = "count")#, color = "red")+
  xlab("Year (Monthly)")+
  ylab("Number of NGS test \n (ThermoFisher)")+
  theme_minimal()

clinical |>
  filter(year_month > "2020-10-01") |>
  group_by(year_month)|>
  summarise(n())|>
  ggplot(aes(x=year_month, y = `n()`)) +
  geom_line()

clinical |> view()

class(clinical$pathological_dx)
summary(clinical$pathological_dx)

clinicalTumorType <- clinical |>
  mutate(diagnosis = as.factor(pathological_dx),
         cancer_type = str_to_lower(tumor_type),
         cancer_type = as.factor(cancer_type),
         cancer_type = case_when(str_detect(cancer_type, "colo|rectal") ~ "Colorectal Cancer",
                                 str_detect(cancer_type, "lung") ~ "Lung Cancer",
                                 str_detect(cancer_type, "glioblastoma|glioma|ependymoma|astrocytoma") ~ "Brain Cancer",
                                 str_detect(cancer_type, "pancrea") ~ "Pancreas Cancer",
                                 str_detect(cancer_type, "gastric") ~ "Stomach Cancer",
                                 str_detect(cancer_type, "breast") ~ "Breast Cancer",
                                 is.na(cancer_type) ~ NA_character_,
                                 TRUE ~ "Other Cancer"),
         cancer_type = as.factor(cancer_type),
         cancer_type = fct_infreq(cancer_type))

clinicalTumorType |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x=year_month)) + #, color = cancer_type) +
  #geom_bar()+
  geom_line(stat = "count")+
  facet_wrap(vars(cancer_type), ncol = 3) + #, scales = "free_y")+
  xlab("Year (Monthly)")+
  ylab("Number of NGS test \n (ThermoFisher)")+
  theme_minimal()

mutation <- readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))

data <- mutation |>
  full_join(clinicalTumorType, by = "pathology_num")

skimr::skim(data)

load(here::here("RData/dfNested.RData"))
dfNested |> str()

dfNested[[2]][[1]]
dfNested[[1]][1]
