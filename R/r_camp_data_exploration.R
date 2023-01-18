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

clinical <- clinical |>
  select(prescription_date, createDate)|>
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
