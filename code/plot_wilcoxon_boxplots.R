## Plot which are different
## Facet by comparison probably
## Maybe make a plot function?

source("code/hypothesis_testing.R")


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

