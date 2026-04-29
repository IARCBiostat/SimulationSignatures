rm(list=ls())
library(ggplot2)
library(tidyverse)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data scenario 1 & 2
res_scenario_1 <- readRDS("results/Res_SecondLayerM_TRUE.rds")
res_scenario_2 <- readRDS("results/res_scenario_2.rds")
source("R/prep_res_scenario_2.R")
# Load data scenario 3
res_scenario_3 <- source("R/prep_res_scenario_3.R")

# Get 
Recap_Ave_s1 <- res_scenario_1 %>% mutate(Corr = replace_na(Corr, 0), Cardinality=CardS, Sensitivity=Sensi, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(cols = c("Cardinality", "Sensitivity", "Specificity", "Correlation"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoM, SecondLayerM, pMrel) %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = "(i)") %>% 
  filter(Method %in% c("lasso.1se"), pMrel == 25, rhoM == 0) %>% 
  dplyr::select(beta, Crit, Mean, Screen, Scenario, pMrel)

Recap_Ave_plot_all <- Recap_Ave_s1 %>% 
  bind_rows(Recap_Ave_s2, Recap_Ave_s3_plots) %>% 
  dplyr::select(beta, Crit, Mean, Screen, Scenario, pMrel)

# Plot for main paper based on Figure 1 scenario i and ii
ggp <- ggplot(Recap_Ave_plot_all, aes(x=beta, y=Mean)) +
  geom_line(aes(linetype=Screen), size=0.7) +
  geom_point() +
  geom_hline(
    data = distinct(Recap_Ave_plot_all, pMrel, Crit, Scenario) %>%
      filter(Crit == "Cardinality") %>%
      mutate(yint = pMrel),
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey"
  ) +
  facet_grid(Crit~Scenario, labeller = label_both, scales = "free_y")+ 
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +  
  labs(
    x = expression(beta),                   
    y = "Average Performance",              
    linetype = "Screening"
  ) +
  theme_bw(base_size = 12) +                 
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(size = 0.3),
    panel.grid.minor = element_blank()
  )
ggsave(ggp, file="results/Figures/scenario_1-3_main.png", width = 9, height = 12, units = "in", dpi=300)
