library(tidyverse)

# - Read expression data
# - Convert from wide to long format
# - Split identifiers to condition, time and replicates

exp_data <- readRDS("data/th_mycn.rds")$exp_by_gene
pheno <- readRDS("data/th_mycn_pheno.rds")

expression_mycn <- exp_data %>%
  pivot_longer(
    -gene,
    names_to = "filename",
    values_to = "expression"
  ) %>%
  mutate(pheno[filename,]) %>%
  select(-filename) %>%
  rename(time = age)

# Make plot
plt_timecourse <- expression_mycn %>%
  # Genes of interest
  filter(gene %in% c( "Ruvbl1", "Ruvbl2")) %>%
  mutate("condition"=fct_recode(genotype, "WT"="TH-MYCN -/-")) %>%
  ggplot(aes(x=time, y=expression, colour=condition, group=condition)) +
  stat_summary(fun = mean, geom = "line") + 
  stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1) +
  facet_wrap(~gene, nrow = 1, scales = "free_y") +
  labs(x = "time (weeks)") +
  ggpubr::theme_pubclean()

plt_timecourse
