rm(list=ls())
library(ggplot2)
library(tidyverse)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data
res_scenario_1 <- readRDS("results/Res_SecondLayerM_FALSE.rds")
res_scenario_2 <- readRDS("results/Res_SecondLayerM_TRUE.rds")

res_all <- rbind(res_scenario_1, res_scenario_2)

# Get 
Recap_Ave <- res_all %>% mutate(Corr = replace_na(Corr, 0), Cardinality=CardS, Sensitivity=Sensi, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(cols = c("Cardinality", "Sensitivity", "Specificity", "Correlation"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoM, SecondLayerM, pMrel) %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = case_when((SecondLayerM==F & rhoM==0) ~ "Baseline",
                              (SecondLayerM==F & rhoM==0.5) ~ "(i)",
                              (SecondLayerM==T & rhoM==0) ~ "(ii)",
                              TRUE ~ "extra"))

count=1
for(slm in c(F,T)){
  for(pm in c(0,0.5)){
    if(count<4){
        ggp <- ggplot(Recap_Ave %>% filter(SecondLayerM == slm, rhoM == pm), aes(x=beta, y=Mean, colour=Method)) +
          geom_line(aes(linetype=Screen), size=0.7) +
          geom_point() +
          facet_grid(Crit~pMrel, labeller = label_both, scales = "free_y")+ 
          # geom_hline(yintercept = pMrel) +
          scale_linetype_manual(values = c("dotted", "solid")) +
          scale_colour_brewer(palette = "Dark2") +  
          labs(
            x = expression(beta),                   
            y = "Average Performance",              
            colour = "Method",                       
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
    # Supplementary figures with all pMrel and 
    ggsave(ggp, file=paste0(paste("results/Figures/AllCrit_Scenario",count,"slm", slm, "metcors", pm, sep = "_"),".png"), width = 9, height = 12, units = "in", dpi=300)}
    count <- count + 1
  }
}

# Plot for main paper based on Figure 1 scenario i and ii
ggp <- ggplot(Recap_Ave %>% filter(pMrel == 25, Scenario != "extra"), aes(x=beta, y=Mean, colour=Method)) +
  geom_line(aes(linetype=Screen), size=0.7) +
  geom_point() +
  facet_grid(Crit~Scenario, labeller = label_both, scales = "free_y")+ 
  # geom_hline(yintercept = pMrel) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +  
  labs(
    x = expression(beta),                   
    y = "Average Performance",              
    colour = "Method",                       
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
ggsave(ggp, file=paste0(paste("results/Figures/AllCrit_pMrel25_main.png"), width = 9, height = 12, units = "in", dpi=300))
