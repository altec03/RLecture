library(tidyverse)

mutation <- readr::read_csv(here::here("data/Pathology_NGS_mutation.csv"))
clinical <- readr::read_csv(here::here("data/Pathology_NGS_1st.csv"), na = "NULL")

library(lubridate)

clinicalDate <- clinical |>
  mutate(age = as.numeric(age),
         prescription_date = ymd(prescription_date),
         year = lubridate::year(prescription_date),
         month = lubridate::month(prescription_date),
         year_month = floor_date(prescription_date, "1 month"))

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

skimr::skim(clinicalDate)

clinicalTumorType <- clinicalDate |>
  mutate(tumor_type = as.factor(tumor_type))

levels(clinicalTumorType$tumor_type)




