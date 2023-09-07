library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

# Store 16S data in datalist-structures (see https://github.com/dermilke/ExCom)

# Project-parameter stores information about location:
# twins-masson = France (various)
# long_term = Frankenhausen, Germany
# mikodu = Darmstadt, Germany
# biodyn_geisenheim = Geisenheim, Germany
# spray = timeseries data (containing France + Frankenhausen)

# Load datalists and add combinatorial parameter for summarizing samples later on
data_16S <- import_data("data/Count_Data_Spain/16S/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  filter_station_datalist(Project %in% c("twins_masson", "long_term", "mikodu", "biodyn_geisenheim", "spray")) %>%
  mutate_meta_datalist(City = ifelse(City == "Kassel", "Frankenhausen", City)) %>%
  mutate_meta_datalist(Kombi2 = ifelse(Country == "France", Country, City)) %>%
  mutate_meta_datalist(Time2 = ifelse(Time == "T0" | is.na(Time), "Early sampling (T0)", "Later sampling (T1)"))

data_ITS <- import_data("data/Count_Data_Spain/ITS/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  filter_station_datalist(Project %in% c("twins_masson", "long_term", "mikodu", "biodyn_geisenheim", "spray")) %>%
  mutate_meta_datalist(City = ifelse(City == "Kassel", "Frankenhausen", City)) %>%
  mutate_meta_datalist(Kombi2 = ifelse(Country == "France", Country, City)) %>%
  mutate_meta_datalist(Time2 = ifelse(Time == "T0" | is.na(Time), "Early sampling (T0)", "Later sampling (T1)"))

# Summarize abundance data by calculating average abundance of ASVs in all
# Location (Kombi2), Category, Time (Time2) - combinations
data_16S_summarized <- data_16S %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  summarize_by_param(Kombi2, Category, Time2)

data_ITS_summarized <- data_ITS %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  summarize_by_param(Kombi2, Category, Time2)

# Transform wide abundance table format into flat table-format by merging ASVs
# on Genus level and group all Genera below 0.9% abundance into "Others"
datatable_16S <- create_datatable(data_16S_summarized, grpBy = Genus, otherThreshold = 0.009,
                                  addColorScheme = T, addOthers = T)

datatable_ITS <- create_datatable(data_ITS_summarized, grpBy = Genus, otherThreshold = 0.009,
                                  addColorScheme = T, addOthers = T)

# Visualize summarized abundance data and taxonomy with stacked barplots
# and change some labels for readability
datatable_16S$table %>%
  mutate(Time2 = ordered(Time2, levels = c("Early sampling (T0)", "Later sampling (T1)"))) %>%
  mutate(Category = ifelse(Category == "Organic+biodynamic", "Biodynamic", Category)) %>%
  mutate(Category = ordered(Category, levels = c("Conventional", "Organic", "Biodynamic"),
                            labels = c("Conv.", "Organic", "Biodynamic"))) %>%
ggplot(., aes(x = Kombi2, y = Abundance * 100, fill = Group)) +
  geom_bar(position = "stack", stat = "identity", col = "black", size = .1) +
  facet_grid(Time2~Category, scales = "free_x", space = "free") +
  scale_fill_manual(values = datatable_16S$color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(.01, "in")) +
  labs(y = "Abundance (%)", x = "", fill = "Class - Genus") + 
  guides(fill = guide_legend(ncol = 1, keywidth = .75, keyheight = .75, byrow = F))

ggsave(filename = "figs/Barplot_Complete_16S.png", width = 8, height = 6.5, dpi = 300)

datatable_ITS$table %>%
  mutate(Time2 = ordered(Time2, levels = c("Early sampling (T0)", "Later sampling (T1)"))) %>%
  mutate(Category = ifelse(Category == "Organic+biodynamic", "Biodynamic", Category)) %>%
  mutate(Category = ordered(Category, levels = c("Conventional", "Organic", "Biodynamic"),
                            labels = c("Conv.", "Organic", "Biodynamic"))) %>%
  ggplot(., aes(x = Kombi2, y = Abundance * 100, fill = Group)) +
  geom_bar(position = "stack", stat = "identity", col = "black", size = .1) +
  facet_grid(Time2~Category, scales = "free_x", space = "free") +
  scale_fill_manual(values = datatable_ITS$color) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(.01, "in")) +
  labs(y = "Abundance (%)", x = "", fill = "Class - Genus") + 
  guides(fill = guide_legend(ncol = 1, keywidth = .75, keyheight = .75, byrow = F))

ggsave(filename = "figs/Barplot_Complete_ITS.png", width = 8, height = 6.5, dpi = 300)
