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

delta_table <- rt_compare %>% 
  select(-PowerMag, -PowerSoil, -Zymobiomics) %>% 
  pivot_longer(cols = c(delta_PMPS, delta_PMZymo, delta_PSZymo),
               names_to = "metric", values_to = "value")



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

