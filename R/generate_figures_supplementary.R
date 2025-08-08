rm(list=ls())
library(ggplot2)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data
Res <- readRDS("results/Res_withLatents_pE_1_SecondLayerL_FALSE_rhoM_0_typeCorrM_allM_sim1_100.rds")

# Prepare dataframe for plotting
Recap_Ave <- Res %>%mutate(Corr = replace_na(Corr, 0), Sensitivity=Sensi, "Latent sensitivity"=Sensi_L, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(cols = c("Sensitivity", "Latent sensitivity", "Specificity", "Correlation"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoL)  %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = case_when(rhoL==0 ~ "Baseline",
                              rhoL==0.5 ~ "Latent scenario"))

# Plot results
ggp <- ggplot(Recap_Ave, aes(x = beta, y = Mean, colour = Method)) +
  geom_line(aes(linetype = Screen), size = 1) +
  geom_point() +
  facet_grid(Crit ~ Scenario, labeller = label_value, scales = "free_y") +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +  
  labs(
    x = expression(beta),                   
    y = "Average Performance",              
    colour = "Method",                       
    linetype = "Screening"
  ) +
  theme_bw(base_size = 14) +                 
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(size = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(plot=ggp, file=paste0("results/Figures/AllCrit_withLatents.png"), width = 10, height = 12, units = "in", dpi=300)