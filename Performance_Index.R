library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")
source("R/get_PGPM.R")

# Store 16S data in datalist-structures (see https://github.com/dermilke/ExCom)

# Project-parameter stores information about location:
# twins-masson = France (various)
# long_term = Frankenhausen, Germany
# mikodu = Darmstadt, Germany
# biodyn_geisenheim = Geisenheim, Germany
# spray = timeseries data (containing France + Frankenhausen)

data_16S <- import_data("data/Count_Data_Spain/16S/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) 

data_ITS <- import_data("data/Count_Data_Spain/ITS/", kingdom = "Prok", abundance_filter = F, min_counts = 2000)

# Add data about hormone-indices to the environmental context data
data_16S$Meta_Data <- data_16S$Meta_Data %>%
  left_join(., read_csv("data/Hormone_Indices.csv"))

data_ITS$Meta_Data <- data_ITS$Meta_Data %>%
  left_join(., read_csv("data/Hormone_Indices.csv"))

# Rarefy datalists and calculate relative abundances
data_16S_mod <- data_16S %>%
  rarefy_datalist(rare_lim = min(.$Meta_Data$Counts_Total)) %>%
  mutate_count_datalist(function(x) x/sum(x))

data_ITS_mod <- data_ITS %>%
  rarefy_datalist(rare_lim = min(.$Meta_Data$Counts_Total)) %>%
  mutate_count_datalist(function(x) x/sum(x))

# Assign plant growth promoting effects to ASVs based on their genus taxonomy
# using in-house curated database (from literature review)
# And retrieve relative abundances of these ASVs and their PGP effects
# We group here into high level "Effects" (Stress Adaptation or Hormone Production)
# and low level "Mechanisms" (detailed PGP Mechanism, such as Auxin production 
# or Phosphate solubilization)
pgpd_list <- get_pgpm(data_16S_mod, data_ITS_mod, richness = F, by_taxonomy = T) 

# Combine functional data with environmental context data and add
# 16S abundance data and ITS abundance data individually
# Also, remove one extreme outlier sample (where total abundance of PGP ASVs was > 37%)
individual_effects <- left_join(select(pgpd_list$Effect_Bac, Sample_ID, Effect, Abundance), 
                                select(pgpd_list$Effect_Fung, Sample_ID, Effect, Abundance), 
                                by = c("Sample_ID", "Effect")) %>%
  dplyr::rename("Abundance_Bacteria" = "Abundance.x", "Abundance_Fungi" = "Abundance.y") %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
  left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
  filter(Project != "Preparation" & Abundance_Bacteria < 0.37)

# Create function that defines linear model that interpolates:
#
# lm(BeCrop Index values ~ Relative abundance of 16S ASVs that perform plant growth promoting effects +
#                          Relative abundance of ITS ASVs that perform plant growth promoting effects)
#
# And then interpolate modelled BeCrop Index values from "measured" relative abundance of PGP ASVs
# Return a table that contains interpolated data, measured data and contextual data together
# with the linear model
test_effect <- function(effect_table, effect, param) {

  model <-  tibble(horm_index = pull(filter(effect_table, Effect == effect), param),
                   horm_measured_Bac = filter(effect_table, Effect == effect)$Abundance_Bacteria*100,
                   horm_measured_Fung = filter(effect_table, Effect == effect)$Abundance_Fungi*100) %>%
    with(., lm(data = ., formula = horm_index ~ horm_measured_Bac + horm_measured_Fung))
  
  table_return <- effect_table %>%
    filter(Effect == effect) %>%
    mutate(Abundance = filter(effect_table, Effect == effect)$Abundance_Bacteria*100 * model$coefficients[2] +
             filter(effect_table, Effect == effect)$Abundance_Fungi*100 * model$coefficients[3])

  return(list(model = model_new, table = table_return))
}

# Apply function to infer model for Hormone induced Biostimulation
model_new <- test_effect(effect_table = individual_effects, effect = "Biostimulation", param = "Hormones")
summary(model_new$model_old)

ggplot(model_new$table, aes(x = Abundance, y = Hormones)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, col = "darkred") +
  theme_bw() +
  labs(y = "BiomeMaker hormone index", x = "Putative hormone producing bacterial and fungi genera") +
  lims(x = c(0.5, 4), y = c(1, 5.5))

ggsave("figs/Hormones_Performance.png", width = 6, height = 5, dpi = 300)

# Apply function to infer model for Stress-Adaptation induced Bioprotection
model_new <- test_effect(effect = "Bioprotection", param = "Stress_Adaptation")
summary(model_new$model_old)

ggplot(model_new$table, aes(x = Abundance, y = Stress_Adaptation)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, col = "darkred") +
  theme_bw() +
  labs(y = "BiomeMaker stress-adaptation index", x = "Putative bioprotecting bacterial and fungi genera") +
  lims(x = c(0.1, 3.5), y = c(1.5, 6))

ggsave("figs/Stress_Performance.png", width = 6, height = 5, dpi = 300)
