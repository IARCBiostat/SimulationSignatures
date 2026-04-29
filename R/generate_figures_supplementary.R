rm(list=ls())
library(ggplot2)
unscreen <- Vectorize(function(char){strsplit(char, "_screen|_noscreen")[[1]][1]})

# Load data ----
## Scenario 1
res_scenario_1 <- readRDS("results/Res_SecondLayerM_TRUE.rds")

# Format scenario 
Recap_Ave <- res_scenario_1 %>% mutate(Corr = replace_na(Corr, 0), Cardinality=CardS, Sensitivity=Sensi, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(cols = c("Cardinality", "Sensitivity", "Specificity", "Correlation"), names_to = "Crit") %>% 
  group_by(Approach, beta, Crit, rhoM, SecondLayerM, pMrel) %>% 
  summarise(Mean=mean(value, na.rm=T),
            SD=sd(value, na.rm=T),
            N=sum(!is.na(value))) %>%          
  mutate(Screen = as_factor(case_when(grepl("noscreen", Approach) ~ "No screen", T ~ "Screen")), 
         Method = unscreen(Approach),
         Scenario = case_when((SecondLayerM==F & rhoM==0.5) ~ "(i)",
                              (SecondLayerM==T & rhoM==0) ~ "(ii)",
                              TRUE ~ "extra")) %>% 
  rename(pChild = pMrel)

# Load scenario 2
source("R/prep_res_scenario_2.R")
# Load scenario 3
source("R/prep_res_scenario_3.R")

# Scenario 1 - full results ----
ggp <- ggplot(Recap_Ave %>% filter(Method == "lasso.1se", SecondLayerM == T, rhoM == 0), aes(x=beta, y=Mean)) +
  geom_line(aes(linetype=Screen), size=0.7) +
  geom_point() +
  geom_hline(
    data = distinct(Recap_Ave, pChild, Crit, Scenario) %>%
      filter(Crit == "Cardinality") %>%
      mutate(yint = pChild),
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey"
  ) +
  facet_grid(Crit~pChild, labeller = label_both, scales = "free_y")+ 
  # geom_hline(yintercept = pChild) +
  scale_linetype_manual(values = c("dotted", "solid")) +
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
ggsave(ggp, file="results/Figures/fig_a1_scen1.png", width = 9, height = 12, units = "in", dpi=300)


# Scenario 2: Supplementary plots ----
# Plot for beta = 0.3; prel = 5,25,125; for "maximum" Bmother, Bstepbrother
## Combinations:
### prel = 5, Bmother/stepbrother = 25/7
### prel = 25, Bmother/stepbrother = 5/7
### prel = 125, Bmother/stepbrother = 3/2
res_plot_2 <- res_plot_scen2 %>% 
  mutate(
    Corr = case_when(is.na(Corr) ~ 0, .default = Corr),
    Comb = paste0("D=", Bhiddenmothers, ", d=", Bstepbrothers)
  ) %>%
  rename(Cardinality=CardS, Sensitivity=Sensi, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(
    cols = c(Cardinality, Sensitivity, Specificity, Correlation, Screen_child, Screen_nonchild),
    names_to = "Crit",
    values_to = "Mean"
  ) %>% 
  mutate(Screen = ifelse(grepl("noscreen", Approach), "No screen", "Screen")) %>%
  filter(
    (prel == 5   & Bhiddenmothers == 25 & Bstepbrothers == 7) |
      (prel == 5   & Bhiddenmothers == 7  & Bstepbrothers == 25) |
      (prel == 25  & Bhiddenmothers == 5  & Bstepbrothers == 7) |
      (prel == 25  & Bhiddenmothers == 7  & Bstepbrothers == 5) |
      (prel == 125 & Bhiddenmothers == 3  & Bstepbrothers == 2) |
      (prel == 125 & Bhiddenmothers == 2  & Bstepbrothers == 3),
    !grepl("Screen", Crit)
  )

ggp2 <- ggplot(res_plot_2 %>% 
                 filter(rho_g == 0.2), 
               aes(x = beta, y = Mean,
                   colour = Comb,
                   linetype = Screen,
                   group = interaction(Comb, Screen))) +
  geom_line(size = 0.7) +
  geom_point() +
  facet_grid(Crit ~ prel, scales = "free_y") +
  geom_hline(
    data = distinct(res_plot_2, prel, Crit) %>%
      filter(Crit == "Cardinality") %>%
      mutate(yint = prel),
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey"
  )+
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = expression(beta),
    y = "Average Performance",
    colour = "D / d",
    linetype = "Screening"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size = 12, face = "bold")
  )
ggp2
ggsave(plot = ggp2, filename = "results/Figures/fig_a3_scen2-rhog0.2.png", height = 12, width = 9)


ggp3 <- ggplot(res_plot_2 %>% 
                 filter(rho_g == 0), 
               aes(x = beta, y = Mean,
                   colour = Comb,
                   linetype = Screen,
                   group = interaction(Comb, Screen))) +
  geom_line(size = 0.7) +
  geom_point() +
  facet_grid(Crit ~ prel, scales = "free_y") +
  geom_hline(
    data = distinct(res_plot_2, prel, Crit) %>%
      filter(Crit == "Cardinality") %>%
      mutate(yint = prel),
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey"
  )+
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = expression(beta),
    y = "Average Performance",
    colour = "D / d",
    linetype = "Screening"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(size = 12, face = "bold")
  )
ggp3
ggsave(plot = ggp3, filename = "results/Figures/fig_a2_scen2-rhog0.png", height = 12, width = 9)

# Plots for the situation prel = 25; varying Bhiddenmothers and Bstepbrothers
plot_metric <- function(df, metric){
  p <- ggplot(df %>% filter(Metric == metric),
              aes(
                x = Bstepbrothers,
                y = factor(Bhiddenmothers),
                fill = Value
              )) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = round(Value, 1)), size = 3, color = "white") +
    facet_grid(Metric ~ Approach, switch = "y") +
    scale_fill_viridis_c(option = "C", na.value = "grey92") +
    theme_minimal()
  return(p)
}

library(patchwork)
res_plot_4 <- res_plot %>% filter(prel == 25, beta==0.3) %>% 
  mutate(Corr = case_when(is.na(Corr) ~ 0, .default = Corr)) %>%
  pivot_longer(
    cols = c(CardS, Sensi, Speci, Corr, Screen_child, Screen_nonchild),
    names_to = "Metric",
    values_to = "Value"
  ) %>% mutate(Screen = ifelse(grepl("noscreen", Approach), "No screen", "Screen")) 

p1 <- plot_metric(res_plot_4 %>% filter(rho_g == 0.2,), "CardS")
p2 <- plot_metric(res_plot_4 %>% filter(rho_g == 0.2,), "Speci")
p3 <- plot_metric(res_plot_4 %>% filter(rho_g == 0.2,), "Sensi")
p4 <- plot_metric(res_plot_4 %>% filter(rho_g == 0.2,), "Corr")

ggp4 <- (p1 + p2) / (p3 + p4)
ggsave(filename = "Figures/nmother_nstep_prel-25_rhog0.2_heatmap.png", plot = ggp4, width=15, height=9)

p1 <- plot_metric(res_plot_4 %>% filter(rho_g == 0), "CardS")
p2 <- plot_metric(res_plot_4 %>% filter(rho_g == 0), "Speci")
p3 <- plot_metric(res_plot_4 %>% filter(rho_g == 0), "Sensi")
p4 <- plot_metric(res_plot_4 %>% filter(rho_g == 0), "Corr")

ggp5 <- (p1 + p2) / (p3 + p4)
ggsave(filename = "Figures/nmother_nstep_prel-25_rhog0_heatmap.png", plot = ggp5, width=15, height=9)

# Scenario 3: Plot all criteria for all n ----

plot_s3_allcrit <- ggplot(Recap_Ave_s3 %>% filter(n %in% c(500, 2500, 62500)) %>% rename(pChild=pMrel), 
                          aes(x=beta, y=value, color = as.factor(n))) +
  geom_line(aes(linetype=Screen), size=0.7) +
  geom_point() +
  facet_grid(
    Crit ~ nfeat,
    scales = "free",
    labeller = labeller(
      nfeat = function(x) paste("p =", x),
      Crit = as_labeller(c(
        "Correlation" = 'bold("Correlation")', 
        "Sensitivity" = 'bold("Sensitivity")',
        "Specificity"     = 'bold("Specificity")',
        "Cardinality"     = 'bold("Cardinality")'
      ), default = label_parsed)
    )
  ) +
  geom_hline(
    data = distinct(Recap_Ave_s3, n, Crit, pMrel, nfeat) %>%
      filter(Crit == "Cardinality") %>% 
      mutate(yint = pMrel),
    aes(yintercept = yint),
    linetype = "dashed",
    colour = "grey"
  ) +
  # geom_hline(yintercept = pChild) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +  
  labs(
    x = expression(beta),                   
    y = "Average Performance",               
    linetype = "Screening",
    color = "Number of samples"
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

ggsave(plot_s3_allcrit, file = "results/Figures/Supp_scen3_allcrit.png", width = 9, height = 12, units = "in", dpi=300)


# Scenario 3: Plot selection frequency for children, descendants, nondescendants ----
all_freq_plot <- all_freq %>% 
  group_by(n, beta, nfeat) %>% 
  summarise(
    child_sel_freq_no_screen  = mean(sel_freq_no_screen[is_child], na.rm = TRUE),
    child_sel_freq_with_screen = mean(sel_freq_with_screen[is_child], na.rm = TRUE),
    desc_sel_freq_no_screen = mean(sel_freq_no_screen[is_desc], na.rm = TRUE),
    desc_sel_freq_with_screen = mean(sel_freq_with_screen[is_desc], na.rm = TRUE),
    nondesc_sel_freq_no_screen = mean(sel_freq_no_screen[is_nondesc], na.rm = TRUE),
    nondesc_sel_freq_with_screen = mean(sel_freq_with_screen[is_nondesc], na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = ends_with(c("no_screen", "with_screen")),
    names_to = c("group", "screening"),
    names_pattern = "(.*)_sel_freq_(.*)",
    values_to = "sel_freq"
  ) %>% mutate(screening=ifelse(screening=="no_screen", "No screen", "Screen"),
               group=ifelse(group=="child", "Child", ifelse(group=="desc", "Descendant", "Non-descendant")))

plot_selfreq <- ggplot(all_freq_plot %>% filter(n %in% c(500, 2500, 62500)), aes(x=beta, y=sel_freq, colour=group)) +
  geom_line(aes(linetype=screening), size=0.7) +
  geom_point() +
  facet_grid(n~nfeat, scales = "free_y",
             labeller = labeller(
               n = function(x) paste("n =", x),
               nfeat = function(x) paste("p =", x))
  ) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_colour_brewer(palette = "Dark2") +  
  labs(
    x = expression(beta),                   
    y = "Average Selection Frequency",              
    colour = "Subset of features",                       
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
ggsave(plot_selfreq, file = "results/Figures/Supp_scen3_selfreq.png", width = 9, height = 12, units = "in", dpi=300)

