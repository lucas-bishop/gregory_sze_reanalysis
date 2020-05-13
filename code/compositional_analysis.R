library(tidyverse)


## basic diversity analysis
shared <- read_tsv("data/mothur_output/final.0.03.subsample.shared", col_types=cols(Group=col_character()))
metadata <- read_csv("data/final.metadata.csv") %>% ifelse(metadata$storage == "ext", mutate())
taxonomy <- read_tsv("data/mothur_output/final.taxonomy") 
  rename_all(tolower) %>% 
  # Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon 
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "_unclassified" = " Unclassified"))) %>% 
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')


otu_data <- shared %>% select(-label, -numOtus) %>%
  pivot_longer(cols=-Group, names_to="otu", values_to="count") %>% 
  mutate(relabund = count / 400)

phylum_data <- inner_join(otu_data, taxonomy)

phylum_taxa <- phylum_data %>%
  group_by(Group, phylum) %>%
  summarize(agg_rel_abund=sum(relabund)) %>% 
  inner_join(., metadata, by = "Group") %>% 
  ungroup()

#pull out top overall genus possibly for strip plot by kit
top_genus <- phylum_data %>% group_by(genus) %>% 
  summarize(agg_rel_abund=sum(relabund)) %>% 
  arrange(desc(agg_rel_abund)) %>% 
  top_n(n=10, agg_rel_abund) %>% 
  pull(genus)






















