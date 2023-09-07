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
data_16S <- import_data("data/Count_Data_Spain/16S/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) 

data_ITS <- import_data("data/Count_Data_Spain/ITS/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) 

# Extract all 16S rRNA gene ASVs that are enriched in BD preparations
# Definition Enriched: if ASV abundance in any preparation > 0.5%
Seq_in_BD_16S <- data_16S %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  filter_station_datalist(!(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD", "Quarzmehl"))) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(apply(select_if(data_16S$Count_Data, is.numeric) %>%
                               select_if(data_16S$Meta_Data$Variety %in% c("500", "500P", "501")) %>%
                               filter(rowSums(.) > 0) %>%
                               mutate_all(function(x) x/sum(x)) %>%
                               mutate_all(function(x) x > 0.005), 1, sum) >= 1) %>%
  .$Count_Data %>%
  .$OTU_ID

# Extract all ITS ASVs that are enriched in BD preparations
# Definition Enriched: if ASV abundance in any preparation > 0.5%
Seq_in_BD_ITS <- data_ITS %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  filter_station_datalist(!(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD", "Quarzmehl"))) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(apply(select_if(data_ITS$Count_Data, is.numeric) %>%
                               select_if(data_ITS$Meta_Data$Variety %in% c("500", "500P", "501")) %>%
                               filter(rowSums(.) > 0) %>%
                               mutate_all(function(x) x/sum(x)) %>%
                               mutate_all(function(x) x > 0.005), 1, sum) >= 1) %>%
  .$Count_Data %>%
  .$OTU_ID

# Calculate relative abundance of count-data and then filter out all ASVs enriched
# in biodynamic preparations
# Then create flat data-table with ASV abundances in all soil-samples
datatable_from_BD_16S <- data_16S %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_16S) %>%
  create_datatable(., grpBy = Genus, otherThreshold = 1,
                   addColorScheme = F, addOthers = T)

datatable_from_BD_ITS <- data_ITS %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_ITS) %>%
  create_datatable(., grpBy = Genus, otherThreshold = 1,
                   addColorScheme = F, addOthers = T)

# Filter out flat data-table for the analysed projects (-> remove samples of pure preparations)
# And filter out samples that have only treated version (-> "Dyn B[1-9]" pattern)
# Then calculate abundance-differences of enriched ASVs between treated and untreated soil samples
datatable_modified_from_BD_16S <- datatable_from_BD_16S %>%
  select(Sample_ID, Abundance) %>%
  left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
  filter(Project %in% c("twins_masson", "biodyn_geisenheim", "spray", "mikodu", "long_term")) %>% 
  filter(!(Project == "long_term" & is.na(Weeks))) %>%
  filter(!str_detect(string = Sample_Reference, pattern = "Dyn B[1-9] - ")) %>%
  select(Block, Project, Sample_ID, Abundance, Weeks, Category, Treatment) %>%
  filter(Category %in% c("Organic+biodynamic", "Organic")) %>%
  group_by(Project, Block, Weeks) %>%
  mutate(Abundance_diff = Abundance[Treatment == "Treatment"] - Abundance[Treatment == "Control"]) %>%
  mutate(Abun_Treat = Abundance[Treatment == "Treatment"]) %>%
  mutate(Abun_Cont = Abundance[Treatment == "Control"]) %>%
  ungroup() %>%
  mutate(Kombi = paste0(Project, "_", Weeks)) %>%
  select(Block, Project, Weeks, Abundance_diff, Abun_Treat, Abun_Cont, Kombi) %>%
  distinct()

datatable_modified_from_BD_ITS <- datatable_from_BD_ITS %>%
  select(Sample_ID, Abundance) %>%
  left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
  filter(Project %in% c("twins_masson", "biodyn_geisenheim", "spray", "mikodu", "long_term")) %>% 
  filter(!(Project == "long_term" & is.na(Weeks))) %>%
  filter(!str_detect(string = Sample_Reference, pattern = "Dyn B[1-9] - ")) %>%
  select(Block, Project, Sample_ID, Abundance, Weeks, Category, Treatment) %>%
  filter(Category %in% c("Organic+biodynamic", "Organic")) %>%
  group_by(Project, Block, Weeks) %>%
  mutate(Abundance_diff = Abundance[Treatment == "Treatment"] - Abundance[Treatment == "Control"]) %>%
  mutate(Abun_Treat = Abundance[Treatment == "Treatment"]) %>%
  mutate(Abun_Cont = Abundance[Treatment == "Control"]) %>%
  ungroup() %>%
  mutate(Kombi = paste0(Project, "_", Weeks)) %>%
  select(Block, Project, Weeks, Abundance_diff, Abun_Treat, Abun_Cont, Kombi) %>%
  distinct()

# Test if abundance of enriched ASVs is significantly higher in treated compared to untreated soils
wilcox.test(datatable_modified_from_BD_16S$Abun_Treat, datatable_modified_from_BD_16S$Abun_Cont,
            alternative = "greater", paired = T)

wilcox.test(datatable_modified_from_BD_ITS$Abun_Treat, datatable_modified_from_BD_ITS$Abun_Cont,
            alternative = "greater", paired = T)

# Plot abundance differences of enriched ASVs between treated and untreated soils
datatable_modified_from_BD_16S %>%
  mutate(Project = ordered(Project, levels = c("biodyn_geisenheim", "twins_masson", "mikodu", "long_term", "spray"),
                           labels = c("Geisenheim\nGemany (n = 4)", "France (n = 21)", "Darmstadt\nGermany (n = 16)",
                                      "Frankenhausen\nGermany (n = 16)", "Time series\nGermany and France (n = 4)"))) %>%
  ggplot(aes(x = as.factor(Weeks), y = Abundance_diff*100)) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 1.2, col = "darkgrey") +
    geom_boxplot() +
    facet_wrap(~Project, scales = "free") +
    theme_bw() +
    labs(y = "Abundance difference of BD prokaryotes\nbetween spray treatment and organic (%)",
         x = "Weeks after spray treatment")

ggsave("figs/BD_organisms_16S_enrichment.png", width = 6, height = 5, dpi = 300)

datatable_modified_from_BD_ITS %>%
  mutate(Project = ordered(Project, levels = c("biodyn_geisenheim", "twins_masson", "mikodu", "long_term", "spray"),
                           labels = c("Geisenheim\nGemany (n = 4)", "France (n = 21)", "Darmstadt\nGermany (n = 16)",
                                      "Frankenhausen\nGermany (n = 16)", "Time series\nGermany and France (n = 4)"))) %>%
  ggplot(aes(x = as.factor(Weeks), y = Abundance_diff*100)) +
    geom_hline(yintercept = 0, linetype = 3, linewidth = 1.2, col = "darkgrey") +
    geom_boxplot() +
    facet_wrap(~Project, scales = "free_x") +
    theme_bw() +
    labs(y = "Abundance difference of fungi\nbetween spray treatment and organic (%)",
         x = "Weeks after spray treatment")

ggsave("figs/BD_organisms_enrichment_ITS.png", width = 6, height = 5, dpi = 300)

#### Barplots of preparation communities ####

# Get microbial composition of biodynamic preparation samples and filter out those taxa
# that were previously identified as "enriched" and transform datalist to
# flat data-table format. 
# Group ASVs by Genus and define "others" as Genera below 1% abundance
datatable_from_BD_16S <- data_16S %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  filter_station_datalist(!(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD", "Quarzmehl"))) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_16S) %>%
  create_datatable(., grpBy = Genus, otherThreshold = 0.01,
                   addColorScheme = T, addOthers = T)

datatable_from_BD_ITS <- data_ITS %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  filter_station_datalist(!(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD", "Quarzmehl"))) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_ITS) %>%
  create_datatable(., grpBy = Genus, otherThreshold = 0.01,
                   addColorScheme = T, addOthers = T)

# Prepare datatable for plotting by changing labels and then
# plot as stacked barplots with colors matching to taxonomy
datatable_from_BD_16S$table %>%
  mutate(City = ifelse(City == "Tobias", "Germany", City)) %>%
  mutate(Variety = paste0("BD", Variety)) %>%
  mutate(Kombi = paste0(City, Variety, Sample_ID)) %>%
  
  ggplot(., aes(x = Kombi, y = Abundance * 100, fill = Group)) +
    geom_bar(position = "stack", stat = "identity", col = "black", size = .1) +
    facet_grid(.~Variety, scales = "free_x", space = "free") +
    scale_fill_manual(values = datatable_from_BD_16S$color) +
    scale_x_discrete(labels = c("Bad Vilbel", "Bad Vilbel", "Zülpich", "Zülpich", "Darmstadt", "Cluny", "Cluny",
                                "Cluny", "Cluny", "Bad Vilbel", "Zülpich", "Darmstadt", "Cluny", "Cluny")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.text = element_text(size = 7.5),
          legend.spacing.y = unit(.01, "in")) +
    labs(y = "Abundance (%)", x = "", fill = "Class - Genus") + 
    guides(fill = guide_legend(ncol = 2, keywidth = .75, keyheight = .75, byrow = F))

ggsave("figs/Barplot_Preparations_16S.png", width = 11, height = 6, dpi = 300)

datatable_from_BD_ITS$table %>%
  mutate(City = ifelse(City == "Tobias", "Germany", City)) %>%
  mutate(Variety = paste0("BD", Variety)) %>%
  mutate(Kombi = paste0(City, Variety, Sample_ID)) %>%
  
  ggplot(., aes(x = Kombi, y = Abundance * 100, fill = Group)) +
    geom_bar(position = "stack", stat = "identity", col = "black", size = .1) +
    facet_grid(.~Variety, scales = "free_x", space = "free") +
    scale_fill_manual(values = datatable_from_BD_ITS$color) +
    scale_x_discrete(labels = c("Bad Vilbel", "Darmstadt", "Cluny", "Cluny",
                                "Cluny", "Cluny", "Bad Vilbel", "Darmstadt", "Cluny", "Cluny")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.text = element_text(size = 7.5),
          legend.spacing.y = unit(.01, "in")) +
    labs(y = "Abundance (%)", x = "", fill = "Class - Genus") + 
    guides(fill = guide_legend(ncol = 1, keywidth = .75, keyheight = .75, byrow = F))

ggsave("figs/Barplot_Preparations_ITS.png", width = 8, height = 6, dpi = 300)

#### NMDS ####

# Calculate Nonmetric multidimensional scaling (NMDS) for count tables of
# biodynamic preparations and their enriched ASVs

# First filter out count-data from biodynamic preparations and next filter out
# previously detected ASVs that are enriched in preparations.
# Next, transform into relative abundance.
datalist_16S_filtered <- data_16S %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_16S)

datalist_ITS_filtered <- data_ITS %>%
  filter_station_datalist(Variety %in% c("Compost_Manure_With_BD", "Compost_Manure_No_BD",
                                         "500", "500P", "501", "Quarzmehl")) %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  filter_taxa_datalist(OTU_ID %in% Seq_in_BD_ITS)

# Transform relative abundance data using "Hellinger" transformation 
# and into Bray-Curtis distances
# prior to calculate nmds using cmdscale()
nmds_16S <- datalist_16S_filtered %>%
  .$Count_Data %>%
  select_if(is.numeric) %>% 
  t() %>%
  vegan::decostand(method = "hellinger") %>%
  vegan::vegdist() %>%
  cmdscale()

nmds_ITS <- datalist_ITS_filtered %>%
  .$Count_Data %>%
  select_if(is.numeric) %>% 
  t() %>%
  vegan::decostand(method = "hellinger") %>%
  vegan::vegdist() %>%
  cmdscale()

# Combine NMDS scores of samples with environmental context-data and visualize
# as scatterplot
datalist_16S_filtered$Meta_Data %>%
  mutate(MDS1 = nmds_16S[,1], MDS2 = nmds_16S[,2]) %>%
  mutate(Variety = ifelse(Variety == "Compost_Manure_No_BD", "Composted Manure\nwithout BD",
                         ifelse(Variety == "Compost_Manure_With_BD", "Composted Manure\nwith BD", Variety)) %>%
           str_replace_all(pattern = "Quarzmehl", "501")) %>%
  
  ggplot(aes(x = MDS1, y = MDS2, fill = Variety, shape = Country)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(21, 22)) +
    scale_fill_manual(values = cbbPalette[c(2, 5, 4, 6, 3)]) +
    labs(fill = "Biodynamic Preparation", shape = "Country") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)),
           shape = guide_legend(override.aes = list(fill = "black", size = 5)),
           size = "none") + 
    theme_bw()

ggsave("figs/NMDS_preparations_16S.png", width = 6, height = 4.5)

datalist_ITS_filtered$Meta_Data %>%
  mutate(MDS1 = nmds_ITS[,1], MDS2 = nmds_ITS[,2]) %>%
  mutate(Variety = ifelse(Variety == "Compost_Manure_No_BD", "Composted Manure\nwithout BD",
                          ifelse(Variety == "Compost_Manure_With_BD", "Composted Manure\nwith BD", Variety)) %>%
           str_replace_all(pattern = "Quarzmehl", "501")) %>%
  
  ggplot(aes(x = MDS1, y = MDS2, fill = Variety, shape = Country)) +
    geom_point(size = 5) +
    scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = cbbPalette[c(2, 5, 4, 6, 3)]) +
    labs(fill = "Biodynamic Preparation", shape = "Country") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)),
           shape = guide_legend(override.aes = list(fill = "black", size = 5)),
           size = "none") +
    theme_bw()

ggsave("figs/NMDS_preparations_ITS.png", width = 6, height = 4.5)
