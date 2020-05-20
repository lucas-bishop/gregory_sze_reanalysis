library(tidyverse)


## basic diversity analysis
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
  mutate(relabund = count / 1000)

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
  top_n(n=12, agg_rel_abund) %>% 
  pull(genus)

genus_taxa <- phylum_data %>% 
  group_by(Group, genus) %>% 
  summarize(agg_rel_abund = sum(relabund)) %>% 
  inner_join(., metadata, by = 'Group') %>% 
  ungroup()

##make strip plot to see if there are any apparent difference between the genera 
top12_stripplot <- genus_taxa %>% filter(genus %in% top_genus) %>% 
ggplot(aes(x= kit, y=agg_rel_abund))+
  geom_hline(yintercept=1/1000, color="gray")+
  geom_jitter(alpha = 0.3, size = 0.2) + 
  geom_boxplot() +
  facet_wrap(~genus) +
  labs(title=NULL, 
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggsave("results/figures/top12_stripchart.png", top12_stripplot, width = 10, height = 6.5)
















