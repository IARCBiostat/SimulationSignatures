BETA <- seq(0.1,0.5,0.1); M17_latent=TRUE

grid <- expand.grid(
  beta = BETA,
  n = c(500, 1000, 2500, 12500, 5*12500, 25*12500),
  M17_latent = M17_latent,
  block = c(55, 110))

## Load data
all_sum  <- tibble()
all_freq <- tibble()

grid <- grid[-c(which(grid$n == 25*12500)),]
for(row in 1:nrow(grid)){
  cat(sprintf("n=%s, beta=%s, block=%s\n", grid[row,"n"], grid[row, "beta"], grid[row, "block"]))
  sim_res <- readRDS(paste0("results/Res_ToyEx_General_n_", grid[row,"n"], "_beta_", grid[row, "beta"], "_M17_latent_TRUE_block_", grid[row, "block"], ".rds"))
  all_sum  <- all_sum  %>% 
    bind_rows(
      as_tibble(sim_res$summary_metrics) %>% 
        mutate(
          n = grid[row,"n"], 
          beta = grid[row, "beta"], 
          nfeat = 19*grid[row, "block"])) 
  all_freq <- all_freq %>% bind_rows(as_tibble(sim_res$freq) %>% mutate(n = grid[row,"n"], beta = grid[row, "beta"], nfeat = 19*grid[row, "block"])) 
}  

all_freq <- all_freq %>% mutate(is_nondesc = case_when(as.numeric(str_remove(feature, "M")) %% 19 == 0 ~ TRUE,
                                                       .default = is_nondesc))

# Prepare dataframe to be compatible w/scenario 1&2 ----
s_s <- c(1,3,6,11,18)
s_nos <- c(1,2,3,6,7,9,11,14,15,18,19)
all_freq_TN <- all_freq %>% 
  group_by(n, beta, nfeat) %>% 
  summarise(
    TP_child_no_screen = 200*mean(sel_freq_no_screen[is_child], na.rm = TRUE),
    TP_child_with_screen = 200*mean(sel_freq_with_screen[is_child], na.rm = TRUE),
    TP_ss_no_screen = 200*mean(sel_freq_no_screen[as.numeric(str_remove(feature, "M"))%%19 %in% s_s], na.rm = TRUE),
    TP_ss_with_screen = 200*mean(sel_freq_with_screen[as.numeric(str_remove(feature, "M"))%%19 %in% s_s], na.rm = TRUE),
    TP_no_screen = 200*mean(sel_freq_no_screen[is_child|is_desc], na.rm = TRUE),
    TP_with_screen = 200*mean(sel_freq_with_screen[is_child|is_desc], na.rm = TRUE),
    TN_no_screen = 200*(1-mean(sel_freq_no_screen[is_nondesc], na.rm = TRUE)),
    TN_with_screen = 200*(1-mean(sel_freq_with_screen[is_nondesc], na.rm = TRUE)),
    FP_no_screen = mean(sel_freq_no_screen[is_nondesc], na.rm = TRUE)*200,
    FP_with_screen = mean(sel_freq_with_screen[is_nondesc], na.rm = TRUE)*200,
    FN_no_screen = 200*(1-mean(sel_freq_no_screen[is_child|is_desc], na.rm = TRUE)),
    FN_with_screen = 200*(1-mean(sel_freq_with_screen[is_child|is_desc], na.rm = TRUE)),
    FN_child_no_screen = 200*(1-mean(sel_freq_no_screen[is_child], na.rm = TRUE)),
    FN_child_with_screen = 200*(1-mean(sel_freq_with_screen[is_child], na.rm = TRUE)),
    FN_ss_no_screen = 200*(1-mean(as.numeric(str_remove(feature, "M"))%%19 %in% s_s, na.rm = TRUE)),
    FN_ss_with_screen = 200*(1-mean(as.numeric(str_remove(feature, "M"))%%19 %in% s_s, na.rm = TRUE))
  ) %>%
  pivot_longer(
    cols = matches("_no_screen$|_with_screen$"),
    names_to = c("metric", "screening"),
    names_pattern = "(.+?)_(no_screen|with_screen)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  )

all_sum_plot <- all_freq_TN %>% left_join(all_sum, by = c("n", "beta", "nfeat", "screening" = "method")) %>%
  mutate(sensi_child = (TP_child)/(TP_child + FN_child), # Sensitivity = TP/(TP+FN)
         speci = (TN)/(TN + FP),
         # sensi_desc  = (TP/(TP + FN)),
         sensi_ss = (TP_ss/(TP_ss + FN_ss))) # Specificity = TN/(TN + FP)

cardi <- all_freq %>% 
  dplyr::select(feature, is_nondesc, nfeat) %>% unique() %>% 
  group_by(nfeat) %>% 
  summarize(n_nondesc = sum(is_nondesc)) %>% 
  mutate(pMrel = case_when(nfeat == 1045 ~ 5*55,
                           nfeat == 2090 ~ 5*110))

Recap_Ave_s3 <- all_sum_plot %>% 
  dplyr::select(n, beta, nfeat, screening, corr.test.noW, n_selected, speci, sensi_child, sensi_ss) %>%
  mutate(Corr = replace_na(corr.test.noW, 0), Cardinality=n_selected,Specificity=speci, Correlation=Corr, `Sensitivity` = sensi_ss) %>% 
  pivot_longer(cols = c("Cardinality", "Specificity", "Correlation", "Sensitivity"), names_to = "Crit") %>% 
  mutate(Screen = as_factor(case_when(grepl("no_screen", screening) ~ "No screen", T ~ "Screen"))) %>% left_join(cardi, by = "nfeat") 

Recap_Ave_s3_plots <- Recap_Ave_s3 %>% filter(nfeat == 1045, n == 1000, !Crit %in% "Sensitivity (S_s)") %>% mutate(Scenario = "(iii)") %>% dplyr::select(beta, Crit, value, Screen, Scenario, pMrel) %>% filter(!grepl("child|desc", Crit)) %>% rename("Mean" = "value")
