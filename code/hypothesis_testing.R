library(tidyverse)
library(broom)
library(purrr)



#hypotesis testing
shared <- read_tsv("data/mothur_output/final.0.03.subsample.shared", col_types=cols(Group=col_character()))
metadata <- read_csv("data/final.metadata.csv")# %>% 
  #filter(kit == "PowerMag" | kit == "PowerSoil" | kit == "Zymobiomics")
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

meta_join <- joined_data %>% full_join(., metadata, by = "Group") %>% drop_na()

# Pick only one storage condition for now - RT


# Need to calculate difference in agg_rel_abund grouped by person, between two kits. With Zero at middle of the x axis

rt_taxa_abund <- meta_join %>% 
  filter(storage == "RoomTemp") %>% 
  group_by(Group, genus) %>% 
  # sum all relative abundances grouped by genera
  mutate(agg_rel_abund = sum(relabund)) %>% 
  select(Group, genus, kit, storage, stool_id, agg_rel_abund) %>% 
  ungroup()

# Made function to create the table with needed p values from hypothesis testing #
# Needs to be paired to be able to look at the change in abundance from one kit to another
wilcoxon_table_genus <- function(data, kit1, kit2){
  data %>% 
  filter(kit == kit1 | kit == kit2) %>% 
    nest(sample_data = c(-genus)) %>%
    mutate(test=map(sample_data, ~tidy(wilcox.test(agg_rel_abund~kit, data=., paired = TRUE)))) %>%
    unnest(test) %>% 
    #mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
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

## Plot which are different
## Facet by comparison probably
## Maybe make a plot function?

PM_PS_plot <- rt_taxa_abund %>% filter(kit == "PowerMag" | kit == "PowerSoil") %>% 
filter(genus %in% sig_PM_PS_genus) %>%
  mutate(genus=factor(genus, levels=sig_PM_PS_genus)) %>%
  mutate(agg_rel_abund=agg_rel_abund+1/21000) %>%
  ggplot(aes(x=genus, y=agg_rel_abund, color=kit)) +
  geom_hline(yintercept=1/1000, color="gray") +
  geom_boxplot(size = 1) +
  geom_point(alpha = 0.1) +
  scale_color_manual(name=NULL,
                      values=c("blue", "red"),
                      breaks=c("PowerMag" , "PowerSoil"),
                      labels=c("PowerMag", "PowerSoil")) +
  labs(x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  theme_classic()

PM_Zymo_plot <- rt_taxa_abund %>% filter(kit == "PowerMag" | kit == "Zymobiomics") %>% 
  filter(genus %in% sig_PM_Zymo_genus) %>%
  mutate(genus=factor(genus, levels=sig_PM_Zymo_genus)) %>%
  mutate(agg_rel_abund=agg_rel_abund+1/21000) %>%
  ggplot(aes(x=genus, y=agg_rel_abund, color=kit)) +
  geom_hline(yintercept=1/1000, color="gray") +
  geom_boxplot(size = 0.8) +
  #geom_point(alpha = 0.1) +
  scale_color_manual(name=NULL,
                     values=c("blue", "green3"),
                     breaks=c("PowerMag" , "Zymobiomics"),
                     labels=c("PowerMag", "Zymobiomics")) +
  labs(x=NULL, y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  coord_flip() +
  theme_classic() + theme(legend.position = "bottom")

# At 1000 subsamples there are no significantly different genera
PS_Zymo_plot <- rt_taxa_abund %>% filter(kit == "PowerSoil" | kit == "Zymobiomics") %>% 
  filter(genus %in% sig_PS_Zymo_genus) %>%
  mutate(genus=factor(genus, levels=sig_PS_Zymo_genus)) %>%
  mutate(agg_rel_abund=agg_rel_abund+1/21000) %>%
  ggplot(aes(x=genus, y=agg_rel_abund, color=kit)) +
  geom_hline(yintercept=1/1000, color="gray") +
  geom_boxplot(size = 1) +
  geom_point(alpha = 0.1) +
  scale_color_manual(name=NULL,
                     values=c("red", "green3"),
                     breaks=c("PowerSoil", "Zymobiomics"),
                     labels=c("PowerSoil", "Zymobiomics")) +
  labs(x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100)) +
  theme_classic()

#myplots <- c(PM_PS_plot, PM_Zymo_plot)

#ggsave("wilccoxon_stripcharts.png", myplots)






