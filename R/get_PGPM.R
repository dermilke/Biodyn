get_pgpm <- function(data_16S, data_ITS, richness = F, by_taxonomy = T) {
  
  pgpb_db <- read_csv("data/database_plant_growth_bacteria.csv") 
  pgpf_db <- read_csv("data/database_plant_growth_fungi.csv")
  
  if (by_taxonomy) {
    merged_bac <- data_16S$Count_Data %>%
      left_join(., select(pgpb_db, Genus, Effect, Mechanism), by = "Genus") %>%
      filter(!is.na(Mechanism))
    
    merged_fung <- data_ITS$Count_Data %>%
      left_join(., select(pgpf_db, Genus, Effect, Mechanism), by = "Genus") %>%
      filter(!is.na(Mechanism))
    
  } else {
    
    merged_bac <- data_16S$Count_Data %>%
      select(-Genus) %>%
      right_join(., select(blasted_16S, -Species), by = "OTU_ID") %>%
      mutate(Genus = Species) %>%
      filter(!is.na(Mechanism))
    
    merged_fung <- data_ITS$Count_Data %>%
      select(-Genus) %>%
      right_join(., select(blasted_ITS, -Species), by = "OTU_ID") %>%
      mutate(Genus = Species) %>%
      filter(!is.na(Mechanism))
    
  }
  
  if (richness) {
    
    merged_bac <- merged_bac %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0))
    
    merged_fung <- merged_fung %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0))
    
    Effect_bac <- merged_bac %>%
      group_by(OTU_ID, Genus, Effect) %>%
      summarize_if(is.numeric, sum) %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0)) %>%
      group_by(Effect) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Effect", "Sample_ID", "Abundance")) %>%
      left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Effect_fung <- merged_fung %>%
      group_by(OTU_ID, Genus, Effect) %>%
      summarize_if(is.numeric, sum) %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0)) %>%
      group_by(Effect) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Effect", "Sample_ID", "Abundance")) %>%
      left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Mechanism_bac <- merged_bac %>%
      group_by(OTU_ID, Genus, Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0)) %>%
      group_by(Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Mechanism", "Sample_ID", "Abundance")) %>%
      left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Mechanism_fung <- merged_fung %>%
      group_by(OTU_ID, Genus, Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0)) %>%
      group_by(Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Mechanism", "Sample_ID", "Abundance")) %>%
      left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
  } else {
    
    Effect_bac <- merged_bac %>%
      group_by(OTU_ID, Genus, Effect) %>%
      summarize_if(is.numeric, mean) %>%
      filter(!is.na(Genus)) %>%
      group_by(Effect) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Effect", "Sample_ID", "Abundance")) %>%
      left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Effect_fung <- merged_fung %>%
      group_by(OTU_ID, Genus, Effect) %>%
      summarize_if(is.numeric, mean) %>%
      filter(!is.na(Genus)) %>%
      group_by(Effect) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Effect", "Sample_ID", "Abundance")) %>%
      left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Mechanism_bac <- merged_bac %>%
      group_by(OTU_ID, Genus, Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      filter(!is.na(Genus)) %>%
      group_by(Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Mechanism", "Sample_ID", "Abundance")) %>%
      left_join(., data_16S$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
    Mechanism_fung <- merged_fung %>%
      group_by(OTU_ID, Genus, Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      filter(!is.na(Genus)) %>%
      group_by(Mechanism) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup() %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("Mechanism", "Sample_ID", "Abundance")) %>%
      left_join(., data_ITS$Meta_Data, by = "Sample_ID") %>%
      as_tibble()
    
  }
  
  return(list(Effect_Bac = Effect_bac, Mechanism_Bac = Mechanism_bac,
              Effect_Fung = Effect_fung, Mechanism_Fung = Mechanism_fung))
}