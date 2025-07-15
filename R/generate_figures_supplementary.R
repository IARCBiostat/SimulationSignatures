rm(list=ls())
library(ggplot2)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data
Res <- readRDS("results/Res_withLatents_pE_1_SecondLayerL_FALSE_rhoM_0_typeCorrM_allM_sim1_100.rds")

Res_nonconverge <- Res %>% 
  group_by(Approach, beta, rhoL) %>% 
  summarise(Count=length(which(is.na(Corr)))) %>%
  mutate(Screening = case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen"),
         Method = unscreen(Approach)) %>%
  pivot_wider(names_from = beta, values_from = Count, names_prefix = "beta=")

# Prepare dataframe for plotting
Recap_Ave <- Res %>% 
  pivot_longer(cols = c("Sensi", "Sensi_L", "Speci", "Corr"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoL)  %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = case_when(rhoL==0 ~ "No shared confounder",
                              rhoL==0.5 ~ "(iii)"))

# Plot results
ggp <- ggplot(Recap_Ave, aes(x = beta, y = Mean, colour = Method)) +
  geom_line(aes(linetype = Screen), size = 1) +
  geom_point() +
  facet_grid(Crit ~ Scenario, labeller = label_both, scales = "free_y") +
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


# ggp_0 <- ggplot(Recap_Ave %>% filter(rhoL==0), aes(x=beta, y=Mean, colour=Method)) +
#   geom_line(aes(linetype=Screen), size=1) +
#   facet_grid(Crit ~ rhoM, labeller = label_both, scales = "free_y")+
#   scale_linetype_manual(values = c("dotted", "solid")) +
#   theme_gray(base_size = 18) +
#   theme(legend.position = "bottom") + ylab("")
# 
# ggsave(plot=ggp_07, file=paste0("Results/Figures/AllCrit_withLatents_rhoL_07.png"), width = 15, height = 15, units = "in", dpi=300)
# ggsave(plot=ggp_0, file=paste0("Results/Figures/AllCrit_withLatents_rhoL_0.png"), width = 15, height = 15, units = "in", dpi=300)
# 
# 
# 
