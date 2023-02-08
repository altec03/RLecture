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

summary(clinicalTumorType$tumor_type)
summary(clinicalDate$tumor_type)

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
         cancer_type = as.factor(cancer_type),
         cancer_type = fct_infreq(cancer_type))

levels(clinicalTumorType$cancer_type)

clinicalTumorType |>
  filter(year_month > "2020-10-01") |>
  ggplot2::ggplot(aes(x=year_month)) + #, color = cancer_type)) +
  #geom_bar()+
  geom_line(stat = "count")+
  facet_wrap(vars(cancer_type))+#, scales = "free")+
  xlab("Year (Monthly)")+
  ylab("Number of NGS test \n (ThermoFisher)")+
  theme_minimal()


mutation <-
  readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))

data <- mutation |>
  full_join(clinicalTumorType, by = "pathology_num")



#######################################################
library(tidyverse)
# load(here::here("RData/data.RData"))
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

geneCount <- data |>
  count(gene)

data <- data |>
  left_join(geneCount)

data <- data |>
  filter(n>20)

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
    values_fill = "wild-type"
  ) |>
  select(-"NA")|>
  left_join(dataTidy, by = "pathology_num")|>
  #colnames()
  pivot_longer(cols = TP53:B2M,
               names_to = "gene",
               values_to = "nucleotide_change") |>
  mutate(nucleotide_change = case_when(nucleotide_change == "wild-type" ~ "wild-type",
                                       nucleotide_change != "wild-type" ~ "mutant"))|>
  group_by(gene)|>
  nest()

dataTidy <- data |>
  select(pathology_num, cancer_type, starts_with("msi"), tumorburden)|>
  distinct_all()


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

#####################################################
library(tidyverse)
load(here::here("RData/dfNested.RData"))

chisq <- dfNested |>
  mutate(msi_chisq = map(data, ~ chisq.test(.x$cancer_type, .x$nucleotide_change)),
         tidy_chisq = map(msi_chisq, broom::tidy))|>
  unnest(tidy_chisq)

chisq |>
  arrange(p.value)|>
  mutate(#gene = as.factor(gene),
         gene = fct_inorder(gene)) |>
  select(gene) |>
  count()

chisq |>
  #arrange(p.value)|>
  #mutate(gene = as.factor(gene),
  #       gene = fct_inorder(gene)) |>
  mutate(significantGene = factor(p.value < 10^-3)) |>
  ggplot(aes(fct_reorder(gene, p.value), -log(p.value, base = 10))) +
  geom_bar(stat = "identity") +
  coord_flip()+
  labs(title = "Differently mutated genes among cancer types",
       subtitle = "Chi-square test",
       x = "",
       y = "Minus log 10 transformed p value")

msi_t.test <- dfNested |>
  unnest()|>
  #ungroup()|>
  mutate(log_msi = log1p(msiscore_number))|>
  group_by(gene)|>
  nest()|>
  mutate(msi_t.test = map(data, ~ t.test(log_msi ~ nucleotide_change,
                                         data = .x)),
         tidy_t.test = map(msi_t.test, broom::tidy)) |>
         #gene = fct_reorder(gene, p.value))|>
  unnest(tidy_t.test)

msi_t.test |>
  mutate(significantGene = factor(p.value < 10^-2)) |>
  ggplot(aes(estimate, -log(p.value, base = 10), label = gene, col = significantGene)) +
  #geom_point() +
  geom_text(size = 2.0)+
  labs(title = "MSI score difference related genes",
       subtitle = "T test",
       x = "Difference",
       y = "Minus log 10 transformed p value")

msi_t.test|>pull(msi_t.test)|>head(2)

colnames(msi_t.test)
dfNested$data[3]
acvr1 <- dfNested$data[dfNested$gene == "ACVR1"][[1]]
acvr1 |>
  mutate(mutation = as.factor(nucleotide_change)) |>
  select(mutation)|>
  summary()

dataMutation|>
  mutate(gene = as.factor(gene))|>
  group_by(cancer_type)|>
  mutate(gene = fct_lump_min(gene, 10))|>
  count(gene)|>
  ungroup()|>
 #mutate(gene = tidytext::reorder_within(gene, n, cancer_type))|>
  filter(gene != "Other")|>
  ggplot(aes(x = gene, y= n))+
  geom_col() +
  #facet_wrap(~cancer_type, scales = "free_y") +
  tidytext::scale_x_reordered() +
  #scale_y_continuous(expand = c(0,0)) +
  coord_flip()

