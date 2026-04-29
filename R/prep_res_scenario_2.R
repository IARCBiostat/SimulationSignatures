source("config/scenario_2.R")
res_scenario_2 <- foreach(rrow=1:nrow(params), .combine="rbind") %do% {
  n = params$n[rrow]
  P = params$P[rrow]
  prel = params$prel[rrow]
  rho_g = params$rho_g[rrow]
  rho_l = params$rho_l[rrow]
  beta = params$beta[rrow]
  Bhiddenmothers = params$Bhiddenmothers[rrow]
  Bstepbrothers = params$Bstepbrothers[rrow]
  iteration=params$iteration[rrow]
  readdir <- "/scratch/wud/simulations/"
  readpath <- paste0(readdir, paste("scenario2",n,P,prel,rho_g,rho_l,beta,Bhiddenmothers,Bstepbrothers,iteration, sep="_"),".rds")
  readRDS(readpath)
}
saveRDS(res_scenario_2, "results/res_scenario_2.rds")

res_plot_scen2 <- res_scenario_2 %>% group_by(Approach, n, P, prel, rho_g, rho_l, beta, Bhiddenmothers, Bstepbrothers) %>%
  summarise(CardS = mean(CardS),
            Sensi = mean(Sensi),
            Speci = mean(Speci),
            Corr = mean(Corr),
            Screen_child = mean(Screen_child),
            Screen_nonchild = mean(Screen_nonchild))

Recap_Ave_s2 <- res_plot_scen2 %>% filter(Bhiddenmothers == 5, Bstepbrothers == 5, prel == 25) %>% 
  mutate(Corr = case_when(is.na(Corr) ~ 0, .default = Corr), Cardinality=CardS, Sensitivity=Sensi, Specificity=Speci, Correlation=Corr) %>%
  pivot_longer(
    cols = c("Cardinality", "Sensitivity", "Specificity", "Correlation"), 
    names_to = "Crit",
    values_to = "Mean"
  ) %>% mutate(Screen = ifelse(grepl("noscreen", Approach), "No screen", "Screen"), Scenario="(ii)") %>%
  filter(rho_g == 0) %>% rename(pMrel=prel) %>%
  dplyr::select(beta, Crit, Mean, Screen, Scenario, pMrel)
