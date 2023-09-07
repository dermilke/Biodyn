library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

# Colorblind-friendly palette
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
color_fun <- colorRampPalette(cbbPalette)

# Store 16S data in datalist-structures (see https://github.com/dermilke/ExCom)

# Project-parameter stores information about location:
# twins-masson = France (various)
# long_term = Frankenhausen, Germany
# mikodu = Darmstadt, Germany
# biodyn_geisenheim = Geisenheim, Germany
# spray = timeseries data (containing France + Frankenhausen)

data_16S <- import_data("data/Count_Data_Spain/16S/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  filter_station_datalist(Project %in% c("twins_masson", "long_term", "mikodu", "biodyn_geisenheim", "spray"))

data_ITS <- import_data("data/Count_Data_Spain/ITS/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  filter_station_datalist(Project %in% c("twins_masson", "long_term", "mikodu", "biodyn_geisenheim", "spray"))

## NMDS
# Transform relative abundance data using "Hellinger" transformation 
# and into Bray-Curtis distances
# prior to calculate nmds using cmdscale()
nmds_16S <- data_16S %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>%
  vegan::decostand(method = "hellinger", MARGIN = 2) %>%
  vegan::vegdist() %>%
  cmdscale()
  
nmds_ITS <- data_ITS %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  t() %>%
  vegan::decostand(method = "hellinger", MARGIN = 2) %>%
  vegan::vegdist() %>%
  cmdscale()

# Prepare NMDS ordination scores for plotting by combining with environmental
# context data and adapt labels for visualization
nmds_16S %>%
  as.tibble() %>%
  magrittr::set_colnames(c("MDS1", "MDS2")) %>%
  mutate(Sample_ID = rownames(nmds_16S)) %>%
  left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
  mutate(City = ifelse(City == "Kassel", "Frankenhausen", City)) %>%
  mutate(City = ifelse(Country == "France", "France - various", City)) %>%
  mutate(Time2 = ifelse(Time == "T0" | is.na(Time), "Early sampling (T0)", "Later sampling (T1)")) %>%
  mutate(Time2 = ordered(Time2, levels = c("Early sampling (T0)", "Later sampling (T1)"))) %>%
  
  ggplot(aes(x = MDS1, y = MDS2, fill = City, shape = Category)) +
    geom_point(size = 2.4) +
    scale_shape_manual(values = c(22, 23, 24, 21)) +
    scale_fill_manual(values = cbbPalette[c(1, 2, 4, 6, 8)]) +
    facet_wrap(~Time2) +
    labs(y = "NMDS2", x = "NMDS2", fill = "Sample location", shape = "Farming practice") +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(size = 3, shape = 21)),
           shape = guide_legend(override.aes = list(size = 3, fill = cbbPalette[5])))

ggsave("figs/NMDS_complete_16S.png", width = 8.5, height = 4.8, dpi = 300)


nmds_ITS %>%
  as.tibble() %>%
  magrittr::set_colnames(c("MDS1", "MDS2")) %>%
  mutate(Sample_ID = rownames(nmds_ITS)) %>%
  left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
  mutate(City = ifelse(City == "Kassel", "Frankenhausen", City)) %>%
  mutate(City = ifelse(Country == "France", "France - various", City)) %>%
  mutate(Time2 = ifelse(Time == "T0" | is.na(Time), "Early sampling (T0)", "Later sampling (T1)")) %>%
  mutate(Time2 = ordered(Time2, levels = c("Early sampling (T0)", "Later sampling (T1)"))) %>%
  
  ggplot(aes(x = MDS1, y = MDS2, fill = City, shape = Category)) +
  geom_point(size = 2.4) +
  scale_shape_manual(values = c(22, 23, 24, 21)) +
  scale_fill_manual(values = cbbPalette[c(1, 2, 4, 6, 8)]) +
  facet_wrap(~Time2) +
  labs(y = "NMDS2", x = "NMDS2", fill = "Sample location", shape = "Farming practice") +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(size = 3, shape = 21)),
         shape = guide_legend(override.aes = list(size = 3, fill = cbbPalette[5])))

ggsave("figs/NMDS_complete_ITS.png", width = 8.5, height = 4.8, dpi = 300)
