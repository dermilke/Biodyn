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
  filter_station_datalist(Project %in% c("twins_masson", "long_term", "mikodu", "biodyn_geisenheim"))

# Add data about hormone-indices to the environmental context data
data_16S$Meta_Data <- data_16S$Meta_Data %>%
  left_join(., read_csv("data/Hormone_Indices.csv"))

# Create table for raw index values 
# and filter out samples from before 1st spray-treatment (long_term, T0) and 
# samples that have not both, treated and untreated variant (Geisenheim, Block 5,6,7,8)
data_red <- data_16S$Meta_Data %>%
  select(Block, Project, Sample_ID, Time, Treatment, 34:45) %>%
  filter(Treatment %in% c("Treatment", "Control")) %>%
  filter(!((Project == "long_term") & (Time == "T0"))) %>%
  filter(!(Project == "biodyn_geisenheim" & Block %in% c("Block_5", "Block_6", "Block_7", "Block_8"))) 

# Transform raw index table from wide data-format to flat data-format
data_red_flat <- map(6:17, function(x) {
    select(data_red, 1:5, x) %>% 
      magrittr::set_names(c(names(data_red[c(1:5)]), "Index")) %>% 
      mutate(Parameter = names(data_red)[x])
  }) %>%
  bind_rows() %>%
  ungroup()

# Test each index within each project if significantly higher in
# treated vs. untreated samples (paired wilcoxon-test)
# and save p-values
test_table <- map(unique(data_red_flat$Project), function(x) {
  
    map(unique(data_red_flat$Parameter), function(y) {
    
      values_index <- data_red_flat %>%
        filter(Project == x, Parameter == y) %>%
        select(Index, Treatment)
      
      test <- wilcox.test(x = pull(filter(values_index, Treatment == "Control"), Index), 
                          y = pull(filter(values_index, Treatment == "Treatment"), Index),
                          paired = TRUE, alternative = "less")
      
      return(tibble(Project = x, Parameter = y, p_val = test$p.value))
    
    }) %>%
    bind_rows()
  }) %>%
  bind_rows() %>%
  mutate(Significant = ifelse(p_val <= 0.05 & p_val > 0.01, "*", 
                              ifelse(p_val <= 0.01 & p_val > 0.001, "**", 
                                     ifelse(p_val < 0.001, "***", 
                                            ifelse(p_val <= 0.1 & p_val > 0.05, "trend", "n.s."))))) %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "_", replacement = "\n") %>%
           str_replace_all(pattern = "\ndiff", "")) %>%
  mutate(Project = ordered(Project, levels = c("long_term", "mikodu", "twins_masson", "biodyn_geisenheim"),
                           labels = c("Frankenhausen\nGermany", "Darmstadt\nGermany", "Various\nFrance", "Geisenheim\nGermany")))

# Plot index-values in boxplots and add significance level
data_red_flat %>%
  mutate(Parameter = str_replace_all(Parameter, pattern = "_", replacement = "\n") %>%
           str_replace_all(pattern = "\ndiff", "")) %>%
  filter(Parameter != "Stress\nSA") %>%
  
  mutate(Project = ordered(Project, levels = c("long_term", "mikodu", "twins_masson", "biodyn_geisenheim"),
                           labels = c("Frankenhausen\nGermany", "Darmstadt\nGermany", "Various\nFrance", "Geisenheim\nGermany"))) %>%
  
  ggplot(aes(x = Treatment, y = Index, fill = Parameter)) +
  geom_boxplot() +
  geom_text(aes(x = "Treatment", y = 5.8, label = Significant), 
            data = filter(test_table, Parameter != "Stress\nSA"),
            size = ifelse(filter(test_table, Parameter != "Stress\nSA")$Significant %in% c("n.s.", "trend"), 3, 7),
            fontface = "bold") +
  scale_fill_manual(values = color_fun(11)) +
  facet_grid(Project~Parameter) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(size = 7),
        legend.position = "none") +
  labs(y = "PGPM Effect Index", x = "", fill = "PGPM Effect")

ggsave("figs/Index_significance.png", width = 9.5, height = 8.5, dpi = 300)

# Prepare dataset for plotting:
# calculate mean-index values for each Project, Parameter, Treatment combination
# and its standard-error
# Then adjust labels for plotting
data_to_plot <- data_red_flat %>%
  mutate(Index = ifelse(Index == 0, NA, Index)) %>%
  group_by(Project, Parameter, Treatment) %>%
  summarize(Index_avg = mean(Index, na.rm = T),
            Index_se = sd(Index, na.rm = T)/sqrt(n())) %>%
  mutate(Treatment = ifelse(Treatment == "Treatment", "Biodynamic", "Control")) %>%
  mutate(Project = str_replace_all(Project, pattern = "long_term", "Frankenhausen") %>%
           str_replace_all(pattern = "biodyn_geisenheim", "Geisenheim") %>%
           str_replace_all(pattern = "mikodu", "Darmstadt") %>%
           str_replace_all(pattern = "twins_masson", "France - various") %>%
           ordered(levels = c("Frankenhausen", "Darmstadt", "France - various", "Geisenheim"))) %>% 
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
  mutate(Treatment = ordered(Treatment, levels = c("Control", "Biodynamic")))

# visualize average index-values and standard-error in barplots
data_to_plot %>%
  ggplot(aes(x = Parameter, y = Index_avg, fill = Treatment)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), col = "black", width = .5, size = .3) +
    geom_errorbar(aes(ymin=ifelse(Index_avg-Index_se < 0, 0, Index_avg-Index_se), ymax=Index_avg+Index_se), width=.2,
                  position=position_dodge(0.8), size = .3) +
    scale_fill_manual(values = c("#548235", "#C00000")) +
    scale_y_continuous(limits = c(0, 5)) +
    scale_x_discrete(expand = c(0.1, 0)) +
    facet_wrap(~Project, nrow = 2) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, size = 7), 
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(y = "BeCrop Index", x = "")

ggsave("figs/Barplot_Index_Diff.png", width = 9, height = 5, dpi = 300)
