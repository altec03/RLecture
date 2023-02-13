library(tidyverse)
load(here::here("RData/dfNested.RData"))

tp53 <- dfNested|>
  pull(data)|>
  purrr::pluck(1)

tp53 <- dfNested|>
  #pull(data)|>
  purrr::pluck("data", 1)

chisq_tp53 <- chisq.test(tp53$cancer_type, tp53$nucleotide_change)

str(chisq_tp53)

chisq_tp53 |> pluck("p.value")

map(dfNested$data, function(x) chisq.test(x$cancer_type, x$nucleotide_change))



chisq <- dfNested |>
  mutate(msi_chisq = map(data, ~ chisq.test(.x$cancer_type, .x$nucleotide_change)),
         p.value = map(msi_chisq, pluck("p.value")))

chisq <- dfNested |>
  mutate(msi_chisq = map(data, ~ chisq.test(.x$cancer_type, .x$nucleotide_change)),
         tidy_chisq = map(msi_chisq, broom::tidy))|>
  unnest(tidy_chisq)

chisq |>
  #mutate(significantGene = factor(p.value < 10^-3)) |>
  ggplot(aes(fct_reorder(gene, p.value), -log(p.value, base = 10))) +
  geom_bar(stat = "identity") +
  coord_flip()+
  theme(axis.text.y=element_text(size = rel(0.6)))+
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

