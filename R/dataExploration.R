library(tidyverse)

clinical <- readr::read_csv(here::here("data/Pathology_NGS_1st.csv"), na = "NULL")

library(lubridate)

clinicalDate <- clinical |>
  mutate(age = as.numeric(age),
         prescription_date = ymd(prescription_date),
         year = lubridate::year(prescription_date),
         month = lubridate::month(prescription_date),
         year_month = floor_date(prescription_date, "bimonth"))

clinicalDate |>
  ggplot(aes(x=age))+
  geom_histogram()

clinicalDate|>
  group_by(year_month)|>
  summarise(n=n(), meanAge = mean(age, na.rm = TRUE))

clinicalDate |>
  ggplot2::ggplot(aes(x=year_month)) +
  geom_bar()

clinicalDate |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x=year_month)) +
  #geom_bar()+
  geom_line(stat = "count")+
  xlab("Year (Monthly)")+
  ylab("Number of NGS test \n (ThermoFisher)")+
  theme_minimal()

clinicalTumorType <- clinicalDate |>
  mutate(tumor_type = as.factor(tumor_type))

levels(clinicalTumorType$tumor_type)

clinicalTumorType <- clinicalDate |>
  mutate(tumor_type = as.factor(tumor_type),
         tumor_type = str_to_lower(tumor_type),
         cancer_type = case_when(str_detect(tumor_type, "colo|rectal") ~ "Colorectal Cancer",
                                 str_detect(tumor_type, "lung") ~ "Lung Cancer",
                                 str_detect(tumor_type, "glioblastoma|glioma|ependymoma|astrocytoma") ~ "Brain Cancer",
                                 str_detect(tumor_type, "pancrea") ~ "Pancreas Cancer",
                                 str_detect(tumor_type, "gastric") ~ "Stomach Cancer",
                                 str_detect(tumor_type, "breast") ~ "Breast Cancer",
                                 is.na(tumor_type) ~ NA_character_,
                                 TRUE ~ "Other Cancer"),
         cancer_type = as.factor(cancer_type))

levels(clinicalTumorType$cancer_type)

clinicalTumorType |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x=year_month)) + #, color = cancer_type)) +
  #geom_bar()+
  geom_line(stat = "count")+
  facet_wrap(vars(cancer_type), ncol = 3)+
  xlab("Year (Monthly)")+
  ylab("Number of NGS test \n (ThermoFisher)")+
  theme_minimal()

mutation <- readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))

dataMutation <- mutation |> 
  left_join(clinicalTumorType, by = "pathology_num")|>
  select(age, gender, pathological_dx, cancer_type, starts_with("tumor"), msiscore, 5:11)

dataMutation|>
  mutate(gene = as.factor(gene))|>
  group_by(cancer_type)|>
  mutate(gene = fct_lump_min(gene, 10))|>
  count(gene)|>
  ungroup()|>
  mutate(gene = tidytext::reorder_within(gene, n, cancer_type))|>
  filter(gene != "Other")|>
  ggplot(aes(x = gene, y= n))+
  geom_col() +
  facet_wrap(~cancer_type, scales = "free_y") +
  tidytext::scale_x_reordered()# +
  #scale_y_continuous(expand = c(0,0)) +
  #coord_flip()
