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
  filter_station_datalist(Project == "spray")

# Add data about hormone-indices to the environmental context data
data_16S$Meta_Data <- data_16S$Meta_Data %>%
  left_join(., read_csv("data/Hormone_Indices.csv"))

# Extract relevant information from environmental context data
data_red <- data_16S$Meta_Data %>%
  select(Block, Project, Sample_ID, Time, Treatment, 34:45) %>%
  filter(Treatment %in% c("Treatment", "Control")) 

# Transform wide data-format to flat data-format
data_red_flat <- map(6:17, function(x) {
    select(data_red, 1:5, x) %>% 
      magrittr::set_names(c(names(data_red[c(1:5)]), "Index")) %>% 
      mutate(Parameter = names(data_red)[x])
  }) %>%
  bind_rows() %>%
  ungroup()

# Prepare dataset for plotting by improving labels
# and calculate index-differences between treated and untreated samples from the same
# Block, Project, Time - combination 
# Then calculate average index-differences and standard-error of index-differences
# for plotting
data_prep <- data_red_flat %>%
  mutate(Treatment = ifelse(Treatment == "Treatment", "Biodynamic", "Control")) %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "Stress_Adaptation", replacement = "Stress\nadaptation") %>%
           str_replace_all(pattern = "^Hormones_", replacement = "") %>%
           str_replace_all(pattern = "Hormones", replacement = "Phytohormone\nproduction") %>%
           str_replace_all(pattern = "^Stress_", replacement = "") %>%
           str_replace_all(pattern = "Exopolys", replacement = "EPS") %>%
           str_replace_all(pattern = "HeavyMetal", replacement = "Heavy metal\nsolubilization") %>%
           str_replace_all(pattern = "Salt", replacement = "Salt tolerance") %>%
           str_replace_all(pattern = "IA", replacement = "Auxin")) %>%
  mutate(Parameter = ordered(Parameter, levels = c("Phytohormone\nproduction", "Auxin", "Cytokinin", "Gibberellin",
                                                   "Stress\nadaptation", "ABA", "ACC", "EPS", "Heavy metal\nsolubilization",
                                                   "SA", "Salt tolerance", "Siderophore"))) %>%
  mutate(Treatment = ordered(Treatment, levels = c("Control", "Biodynamic")))  %>%
  group_by(Parameter, Time, Block) %>%
  mutate(Index = Index[Treatment == "Biodynamic"] - Index[Treatment == "Control"]) %>%
  filter(Treatment == "Biodynamic") %>%
  group_by(Time, Parameter) %>%
  summarize(Index_avg = mean(Index, na.rm = T),
            Index_se = sd(Index, na.rm = T)/sqrt(n())) %>%
  mutate(Time = ordered(Time, levels = c("T0", "T1", "T2", "T3", "T4", "T5", "T6"),
                        labels = c("0", "2", "4", "6", "8", "11", "15")))

# Plot average index-differences of timeseries in barplots
data_prep %>%   
  ggplot(aes(x = Time, y = Index_avg, fill = Parameter)) +
    geom_bar(stat = "identity", position = "dodge", col = "black", size = .3) +
    geom_linerange(aes(ymin=Index_avg, ymax=ifelse(Index_avg < 0, Index_avg-Index_se, Index_avg+Index_se)), size = .3) +
    facet_wrap(~Parameter, nrow = 2) +
    scale_fill_manual(values = color_fun(20)[c(1, 3, 4, 6, 9, 10, 12, 14, 15, 16, 18, 20)]) +
    labs(y = "BeCrop Index\nDifference between biodynamic and control", x = "Weeks after first spray",
         fill = "") +
    theme_bw() +
    theme(legend.position = "none")

ggsave("figs/Barplot_Timeseries_complete_errorbar.png", width = 9, height = 5, dpi = 300)

# Plot average across all index-differences of timeseries in barplots
data_prep %>%
  group_by(Time) %>%
  summarize(Index = mean(Index_avg),
            Index_se = sd(Index_avg, na.rm = T)/sqrt(n())) %>%
  mutate(Parameter = "Average of all soil health indices") %>%
  ggplot(aes(x = Time, y = Index)) +
  geom_bar(stat = "identity", position = "dodge", col = "black", size = .3, fill = cbbPalette[4]) +
  geom_linerange(aes(ymin=Index, ymax=ifelse(Index < 0, Index-Index_se, Index+Index_se)), size = .3) +
  labs(y = "BeCrop Index\nDifference between biodynamic and control", x = "Weeks after first spray") +
  theme_bw() +
  facet_wrap(~Parameter)

ggsave("figs/Barplot_Timeseries_avg_errorbar.png", width = 4, height = 4, dpi = 300)
