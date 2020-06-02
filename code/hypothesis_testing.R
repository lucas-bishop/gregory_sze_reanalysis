library(tidyverse)
library(broom)
library(purrr)



#hypotesis testing
shared <- read_tsv("data/mothur_output/final.0.03.subsample.shared", col_types=cols(Group=col_character()))
# Just want one storage condition for now
metadata <- read_csv("data/final.metadata.csv") %>% #filter(storage == "RoomTemp") %>% 
  filter(kit == "PowerMag" | kit == "PowerSoil" | kit == "Zymobiomics")
taxonomy <- read_tsv("data/mothur_output/final.taxonomy") %>% 
  rename_all(tolower) %>% 
  # Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon 
                                              'Bacteria_unclassified' = 'Unclass.',
                                              "_unclassified" = " Unclass."))) %>% 
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')


otu_data <- shared %>% select(-label, -numOtus) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count") %>% 
  mutate(relabund = count / 1000) 


joined_data <- inner_join(otu_data, taxonomy)


# Need to calculate difference in agg_rel_abund grouped by person, between two kits. With Zero at middle of the x axis

rt_taxa_abund <- joined_data %>%
  group_by(Group, genus) %>%
  # sum all relative abundances grouped by genera
  summarize(agg_rel_abund = sum(relabund)) %>% 
  #filter(agg_rel_abund > 0.000) %>% 
  # need to get rid of rows that contain genera that are different between two kits on the same stool_id
  # need to add a mutate(delta_abund = )
  inner_join(., metadata, by = "Group") %>%
  ungroup() %>% 
  mutate(agg_rel_abund = ifelse(agg_rel_abund == 0, NA, agg_rel_abund))

relevant_genera <- rt_taxa_abund %>% group_by(genus) %>% 
  filter(agg_rel_abund > 0.000) %>% count() %>% 
  mutate(frequency = n/276) %>% filter(frequency >= 0.25) %>% pull(genus)


# convert to wide data frame to set up for calculating deltas
rt_compare <- rt_taxa_abund %>% 
  filter(storage == "RoomTemp", genus %in% relevant_genera) %>% 
  select(-Group, -storage) %>% 
  pivot_wider(id_cols = c(genus, stool_id), names_from = kit, values_from = agg_rel_abund) %>%  
  #drop the rows for stool samples that may have been dropped from one kit
  group_by(genus) %>% 
  mutate(delta_PMPS = ifelse((PowerMag - PowerSoil) == 0, NA, PowerMag - PowerSoil),
         delta_PMZymo = ifelse((PowerMag - Zymobiomics) == 0, NA, PowerMag - Zymobiomics),
         delta_PSZymo = ifelse((PowerSoil - Zymobiomics) == 0, NA, PowerSoil - Zymobiomics))

delta_table <- rt_compare %>% 
  select(-PowerMag, -PowerSoil, -Zymobiomics) %>% 
  pivot_longer(cols = c(delta_PMPS, delta_PMZymo, delta_PSZymo),
               names_to = "metric", values_to = "value")

# get significantly different taxa
wilcoxon_table_genus <- function(data, kit1, kit2){
  data %>% 
  filter(storage == "RoomTemp", genus %in% relevant_genera,kit == kit1 | kit == kit2) %>% 
  select(-Group, -storage) %>%
  pivot_wider(id_cols = c(genus, stool_id), names_from = kit, values_from = agg_rel_abund) %>% 
  drop_na() %>% pivot_longer(cols=-c(genus, stool_id),names_to = "kit", values_to = "agg_rel_abund") %>% 
  nest(sample_data = c(-genus)) %>%
  mutate(test=map(sample_data, ~tidy(wilcox.test(agg_rel_abund~kit, data=., paired = TRUE, alternative = "two.sided")))) %>%
  unnest(test) %>% mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj)
}






# PM - PS comparison
PM_PS_tests <- wilcoxon_table_genus(rt_taxa_abund, "PowerMag", "PowerSoil")

sig_PM_PS_genus <- PM_PS_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)


# PM - Zymo comparison
PM_Zymo_tests <- wilcoxon_table_genus(rt_taxa_abund, "PowerMag", "Zymobiomics")

sig_PM_Zymo_genus <- PM_Zymo_tests %>% 
  filter(p.value.adj <= 0.05) %>%
  pull(genus)

# PS - Zymo comparison
PS_Zymo_tests <- wilcoxon_table_genus(rt_taxa_abund, "PowerSoil", "Zymobiomics")

sig_PS_Zymo_genus <- PS_Zymo_tests %>% 
  filter(p.value.adj <= 0.05) %>% 
  pull(genus)






