rm(list=ls())
library(ggplot2)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data
res_scenario_1 <- readRDS("results/Res_SecondLayerM_FALSE.rds")
res_scenario_2 <- readRDS("results/Res_SecondLayerM_TRUE.rds")

res_all <- rbind(res_scenario_1, res_scenario_2)

# Get 
Recap_Ave <- res_all %>% 
    pivot_longer(cols = c("CardS", "Sensi", "Speci", "Corr"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoM, SecondLayerM, pMrel) %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = case_when((SecondLayerM==F & rhoM==0) ~ "(i)",
                              (SecondLayerM==F & rhoM==0.5) ~ "(ii)",
                              (SecondLayerM==T & rhoM==0) ~ "(iii)",
                              TRUE ~ "extra"))

Recap_Diff <- Recap_Ave %>%
  ungroup() %>%
  select(-Approach) %>%
  pivot_wider(names_from = Screen, values_from = c(Mean, SD, N)) %>%
  group_by(beta, Crit, rhoM, SecondLayerM, pMrel, Method, Scenario) %>%
  summarise(
    se_diff = sqrt(`SD_Screen`^2 / `N_Screen` + `SD_No screen`^2 / `N_No screen`),
    mean_diff = `Mean_Screen` - `Mean_No screen`,
    z = mean_diff / se_diff,
    p = 2 * pnorm(-abs(z)),
    Mean_Screen = `Mean_Screen`,
    .groups = "drop"
  ) %>%
  mutate(Signif = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  ))


Res_nonconverge <- res_all %>% 
  group_by(Approach, beta, rhoM, SecondLayerM, pMrel) %>% 
  summarise(Count=length(which(is.na(Corr)))) %>%
  mutate(Screening = case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen"),
         Method = unscreen(Approach)) %>%
  pivot_wider(names_from = beta, values_from = Count, names_prefix = "beta=")

count=1
for(slm in c(F,T)){
  for(pm in c(0,0.5)){
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
    ggsave(ggp, file=paste0(paste("results/Figures/AllCrit_Scenario",count,"slm", slm, "metcors", pm, sep = "_"),".png"), width = 9, height = 12, units = "in", dpi=300)
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

# Add significance labels with color by l1se/lmin group
## Filter Recap_Diff for pMrel == 25 and Scenario in main plot (if Scenario exists)
sig_labels <- Recap_Diff %>%
  filter(pMrel == 25, Scenario != "extra") %>%
  mutate(Signif = ifelse(Signif == "", NA, Signif))  # Only keep labels where significant

ggp + 
  geom_text(
    data = sig_labels, 
    aes(
      x = beta, 
      y = Mean_Screen,   # position near difference mean
      label = Signif,
      colour = Method
    ), 
    size = 5,
    vjust = -0.5
  ) +
  scale_colour_manual(
    name = "Significance Group",
    values = c("lasso.1se" = "#1B9E77FF", "lasso.min" = "#D95F02FF", "Other" = "grey")
  )


ggsave(ggp, file=paste0(paste("results/Figures/AllCrit_All_Scenario_
                              metcors", pm, sep = "_"),".png"), width = 9, height = 12, units = "in", dpi=300)

