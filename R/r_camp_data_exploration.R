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
  select(prescription_date, createDate)|>
  mutate(prescription_date = lubridate::ymd(prescription_date)) |>
  mutate(year = year(prescription_date))

clinical |>
  select(prescription_date, createDate)|>
  mutate(prescription_date = lubridate::ymd(prescription_date)) |>
  mutate(presc_year = year(prescription_date),
         presc_month = month(prescription_date),
         presc_week = week(prescription_date))|>
  view()
