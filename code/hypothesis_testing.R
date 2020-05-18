library(tidyverse)
library(broom)
library(purrr)



#hypotesis testing
shared <- read_tsv("data/mothur_output/final.0.03.subsample.shared", col_types=cols(Group=col_character()))
metadata <- read_csv("data/final.metadata.csv") %>% 
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
  mutate(relabund = count / 400)

joined_data <- inner_join(otu_data, taxonomy)

#hypothesis testing to see if any are different across the kits or storage conditions

top_genus <- joined_data %>% group_by(genus) %>% 
  summarize(agg_rel_abund=sum(relabund)) %>% 
  arrange(desc(agg_rel_abund)) %>% 
  top_n(n=9, agg_rel_abund) %>% 
  pull(genus)

genus_taxa <- joined_data %>% 
  group_by(Group, genus) %>% 
  summarize(agg_rel_abund = sum(relabund)) %>% 
  inner_join(., metadata, by = 'Group') %>% 
  ungroup()

genus_tests <- genus_taxa %>% nest(sample_data = c(-genus)) %>%
  mutate(test=map(sample_data, ~tidy(wilcox.test(agg_rel_abund~kit, data=.)))) %>%
  unnest(test)

















