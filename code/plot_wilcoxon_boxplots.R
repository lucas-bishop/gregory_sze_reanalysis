## Plot which are different
## Facet by comparison probably
## Maybe make a plot function?

source("code/hypothesis_testing.R")

# convert to wide data frame to set up for calculating deltas

# recalc with change in abundance divided by starting abundance: (PM -PS) / PM
rt_compare <- rt_taxa_abund %>% 
  filter(storage == "RoomTemp", genus %in% relevant_genera) %>% 
  select(-Group, -storage) %>% 
  pivot_wider(id_cols = c(genus, stool_id), names_from = kit, values_from = agg_rel_abund) %>%  
  group_by(genus) %>% 
  mutate(delta_PMPS = ifelse((PowerMag - PowerSoil) == 0, NA, PowerMag - PowerSoil),
         delta_PMZymo = ifelse((PowerMag - Zymobiomics) == 0, NA, PowerMag - Zymobiomics),
         delta_PSZymo = ifelse((PowerSoil - Zymobiomics) == 0, NA, PowerSoil - Zymobiomics))


## Plotting

### Need to look at how many reads were even in these, and gram negative vs gram postive


PM_PS_plot <- rt_compare %>%
  select(genus, stool_id, delta_PMPS) %>% drop_na() %>% 
  filter(genus %in% sig_PM_PS_genus) %>%
  ggplot(aes(x=genus, y=delta_PMPS*100)) +
  geom_boxplot(size = 1) +
  labs(x=NULL,
       y="Change in relative abundance (%)\nfrom PowerMag to PowerSoil") +
  coord_flip() +
  theme_classic()

PM_Zymo_plot <- rt_compare %>% 
  select(genus, stool_id, delta_PMZymo) %>% drop_na() %>% 
  filter(genus %in% sig_PM_Zymo_genus) %>%
  ggplot(aes(x=genus, y=delta_PMZymo*100)) +
  geom_boxplot(size = 1) +
  labs(x=NULL,
       y="Change in relative abundance (%)\nfrom PowerMag to Zymobiomics") +
  coord_flip() +
  theme_classic()

PS_Zymo_plot <- rt_compare %>% 
  select(genus, stool_id, delta_PSZymo) %>% drop_na() %>% 
  filter(genus %in% sig_PS_Zymo_genus) %>%
  ggplot(aes(x=genus, y=delta_PSZymo*100)) +
  geom_boxplot(size = 1) +
  labs(x=NULL,
       y="Change in relative abundance (%)\nfrom PowerSoil to Zymobiomics") +
  coord_flip() +
  theme_classic()



ggsave("delta_PMPS_boxplot.png", PM_PS_plot)
ggsave("delta_PMZymo_boxplot.png", PM_Zymo_plot)
ggsave("delta_PSZymo_boxplot.png", PS_Zymo_plot)

