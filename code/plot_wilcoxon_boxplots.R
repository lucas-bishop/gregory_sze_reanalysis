## Plot which are different
## Facet by comparison probably
## Maybe make a plot function?

source("code/hypothesis_testing.R")

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

# create long version of df for making combined figure
delta_table <- rt_compare %>% 
  select(-PowerMag, -PowerSoil, -Zymobiomics) %>% 
  pivot_longer(cols = c(delta_PMPS, delta_PMZymo, delta_PSZymo),
               names_to = "metric", values_to = "value") %>% 
  drop_na()


## Plotting


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

