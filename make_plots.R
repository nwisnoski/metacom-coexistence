library("tidyverse")

theme_set(theme_light() + 
            theme(strip.text = element_text(color = "black", size = 14),
                  strip.background = element_rect(fill = "gray90"),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  axis.ticks.length = unit(x = .2, units = "cm"),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 12)))

# competition scenarios
#mixed <- "out/last_t_above_total_dyn_mixed_2020-10-05_105308.csv"
equal_above <- "sim_output/final_equal_above.csv"
stable_above <- "sim_output/final_stable_above.csv"
equal_below <- "sim_output/final_equal_below.csv"
stable_below <- "sim_output/final_stable_below.csv"

# below
#mixed <- "out/last_t_below_total_dyn_mixed_2020-09-22_123927.csv"
#equal <- "out/last_t_below_total_dyn_equal_2020-09-22_123617.csv"
#stable <- "out/last_t_below_total_dyn_stabilizing_2020-09-22_124640.csv"

file_list <- c(equal_above, stable_above, equal_below, stable_below)
competition = file_list[1]
make_plots <- function(competition){
  condition_type <- str_remove(competition, ".csv") %>% 
    str_remove("sim_output/final_") %>% 
    str_split("_") %>% .[[1]] %>% .[1]
  
  position <- str_remove(competition, ".csv") %>% 
    str_remove("sim_output/final_") %>% 
    str_split("_") %>% .[[1]] %>% .[2]
  
  
  # read data and ignore s=0.1 condition
  last_t_out <- read_csv(competition) %>% 
    filter(survival != 0.1)
  
  
  disp_rates <- sort(unique(last_t_out$dispersal))
  germ_fracs <- sort(unique(last_t_out$germination))
  surv_fracs <- sort(unique(last_t_out$survival))
  
  alpha_out <- last_t_out %>%
    group_by(dispersal, germination, survival, patch, rep) %>% 
    mutate(N_pa = 1*(N > 0)) %>%
    summarize(alpha = sum(N_pa)) %>%
    ungroup()
  
  alpha_out <- alpha_out %>%
    group_by(dispersal, germination, survival, rep) %>%
    summarize(mean_alpha = mean(alpha))
  
  gamma_out <- last_t_out %>%
    ungroup() %>%
    mutate(N_pa = 1*(N > 0)) %>%
    group_by(dispersal, germination, survival, rep) %>%
    summarize(gamma = length(unique(species)))
  
  div_part <- alpha_out %>%
    left_join(gamma_out) %>%
    mutate(beta = gamma / mean_alpha)
  
  no_sb <- div_part %>% 
    filter(germination == 1) %>% 
    mutate(survival = "No seed bank")
  
  div_part <- div_part %>% 
    mutate(
      survival = case_when(
        survival == 0.5 ~ "Surv = 0.5",
        survival == 1 ~ "Surv. = 1"
      )
    ) %>% bind_rows(no_sb)
  
  div_part %>%
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"),
                              labels = c("mean alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = partition, fill = partition)) +
    geom_point(alpha = 0.05) +
    #geom_line(mapping = aes(group = paste(partition, rep)), stat = "smooth", alpha = 0.2, se = F) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(germination), cols = vars(survival), drop = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "top", axis.text = element_text(size = 8)) +
    ggsave(paste0("figures/diversity_partitioning_", position, "_", condition_type,".pdf"), width = 8, height = 12)
  
  div_part %>%
    #filter(germination == 1) %>% 
    # mutate(#germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2))),
    #   survival = factor(survival, levels = surv_fracs, labels = c("No seed bank", "Surv. = 0.5", "Surv. = 1"))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"),
                              labels = c("mean alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = germination, fill = germination, group = germination)) +
    geom_point(alpha = 0.05) +
    #geom_line(mapping = aes(group = paste(germination, rep)), stat = "smooth", alpha = 0.2, se = F) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(partition), cols = vars(survival), drop = FALSE, scales = "free_y") +
    scale_x_log10() +
    #scale_y_log10() +
    scale_color_viridis_c("Germination", option = "B", end = .8) +
    scale_fill_viridis_c("Germination", option = "B", end = .8) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right") +
    ggsave(paste0("figures/germination_", position, "_", condition_type,".pdf"), width = 12, height = 8)
  
  div_part %>%
    filter(germination <= 0.5) %>% 
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"),
                              labels = c("mean alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = survival, fill = survival, group = survival)) +
    geom_point(alpha = 0.05) +
    #geom_line(mapping = aes(group = paste(survival, rep)), stat = "smooth", alpha = 0.2, se = F) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(partition), cols = vars(germination), drop = TRUE, scales = "free_y") +
    scale_x_log10() +
    #scale_y_log10() +
    scale_color_viridis_d("Survival", option = "D", end = .9) +
    scale_fill_viridis_d("Survival", option = "D", end = .9) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right") +
    ggsave(paste0("figures/survival_", position, "_", condition_type,".pdf"), width = 12, height = 8)
  
  germ_comp <- div_part %>% 
    filter(germination == 1 | germination == 0.1) %>% 
    filter(survival == 1) %>% 
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2))),
           survival = factor(survival, levels = surv_fracs, labels = paste("Surv. =", round(surv_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"), 
                              labels = c("mean alpha", "beta", "gamma"))) 
  # germ_comp %>% 
  #   ggplot(aes(x = dispersal, y = diversity, color = germination, fill = germination)) +
  #   geom_point(alpha = 0.05) +
  #   #geom_line(mapping = aes(group = paste(germination, rep)), stat = "smooth", alpha = 0.2, se = F) +
  #   geom_smooth(alpha = 0.25) +
  #   facet_grid(rows = vars(partition), drop = T, scales = "free_y") +
  #   scale_x_log10() +
  #   #scale_y_log10() +
  #   scale_color_viridis_d("Germination", option = "D", end = .8, begin = .2) +
  #   scale_fill_viridis_d("Germination", option = "D", end = .8, begin = .2) +
  #   labs(y = "Diversity", x = "Dispersal") +
  #   theme(legend.position = "top") +
  #   ggsave(paste0("figures/single_diversity_partitioning_", position, "_", condition_type,".pdf"), height = 6, width = 4)
  # 
  # germ_comp %>% 
  #   ggplot(aes(x = dispersal, y = diversity, color = germination, fill = germination)) +
  #   geom_point(data = filter(germ_comp, germination == "Germ. = 1"), alpha = 0.05) +
  #   #geom_line(mapping = aes(group = paste(germination, rep)), stat = "smooth", alpha = 0.2, se = F) +
  #   geom_smooth(data = filter(germ_comp, germination == "Germ. = 1"), alpha = 0.25) +
  #   facet_grid(rows = vars(partition), drop = T, scales = "free_y") +
  #   scale_x_log10() +
  #   #scale_y_log10() +
  #   scale_color_viridis_d("Germination", option = "D", end = .8, begin = .8) +
  #   scale_fill_viridis_d("Germination", option = "D", end = .8, begin = .8) +
  #   labs(y = "Diversity", x = "Dispersal") +
  #   theme(legend.position = "top") +
  #   ggsave(paste0("figures/single_diversity_partitioning_nodorm_", position, "_", condition_type,".pdf"), height = 6, width = 4)
  
}

### Loop through files
for(competition in file_list){
  make_plots(competition)
}



# look at population abundances
above_stable <- read_csv(stable_above)
above_stable %>% 
  group_by(dispersal, species) %>% 
  summarize(mean_N = mean(N), max_N = max(N), min_N = min(N), sd_N = sd(N)) %>% 
  ggplot(aes(y = mean_N, x = dispersal)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_ribbon(aes(ymax = max_N, ymin = min_N), alpha = 0.2) + 
  scale_x_log10() +
  facet_wrap(~species) +
  theme(axis.text = element_text(size = 8)) + 
  labs(x = "Dispersal rate", y = "Mean abundance\n(+/- max/min)", title = "Aboveground, stable competition") + 
  ggsave("figures/mean_abund_species_above_stable.pdf", width = 10, height = 3/4*10)

above_equal <- read_csv(equal_above)
above_equal %>% 
  group_by(dispersal, species) %>% 
  summarize(mean_N = mean(N), max_N = max(N), min_N = min(N), sd_N = sd(N)) %>% 
  ggplot(aes(y = mean_N, x = dispersal)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_ribbon(aes(ymax = max_N, ymin = min_N), alpha = 0.2) + 
  #geom_ribbon(aes(ymax = mean_N + sd_N, ymin = mean_N - sd_N), alpha = 0.2) + 
  scale_x_log10() +
  facet_wrap(~species) +
  theme(axis.text = element_text(size = 8)) +
  labs(x = "Dispersal rate", y = "Mean abundance\n(+/- max/min)", title = "Aboveground, equal competition") +
  ggsave("figures/mean_abund_species_above_equal.pdf", width = 10, height = 3/4*10)


below_stable <- read_csv(stable_below)
below_stable %>% 
  group_by(dispersal, species) %>% 
  summarize(mean_N = mean(D), max_N = max(D), min_N = min(D), sd_N = sd(D)) %>% 
  ggplot(aes(y = mean_N, x = dispersal)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_ribbon(aes(ymax = max_N, ymin = min_N), alpha = 0.2) + 
  scale_x_log10() +
  facet_wrap(~species) +
  theme(axis.text = element_text(size = 8)) + 
  labs(x = "Dispersal rate", y = "Mean belowground abundance\n(+/- max/min)", title = "Belowground, stable competition") + 
  ggsave("figures/mean_abund_species_below_stable.pdf", width = 10, height = 3/4*10)

below_equal <- read_csv(equal_below)
below_equal %>% 
  group_by(dispersal, species) %>% 
  summarize(mean_N = mean(D), max_N = max(D), min_N = min(D), sd_N = sd(D)) %>% 
  ggplot(aes(y = mean_N, x = dispersal)) + 
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_ribbon(aes(ymax = max_N, ymin = min_N), alpha = 0.2) + 
  #geom_ribbon(aes(ymax = mean_N + sd_N, ymin = mean_N - sd_N), alpha = 0.2) + 
  scale_x_log10() +
  facet_wrap(~species) +
  theme(axis.text = element_text(size = 8)) +
  labs(x = "Dispersal rate", y = "Mean belowground abundance\n(+/- max/min)", title = "Belowground, equal competition") +
  ggsave("figures/mean_abund_species_below_equal.pdf", width = 10, height = 3/4*10)


