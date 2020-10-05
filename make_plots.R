library("tidyverse")
library("gganimate")

theme_set(theme_light() + 
            theme(strip.text = element_text(color = "black"),
                  strip.background = element_rect(fill = "gray90")))

# competition scenarios
mixed <- "out/last_t_above_total_dyn_mixed_2020-09-22_123927.csv"
equal <- "out/last_t_above_total_dyn_equal_2020-09-22_123617.csv"
stable <- "out/last_t_above_total_dyn_stabilizing_2020-09-22_150857.csv"

# below
#mixed <- "out/last_t_below_total_dyn_mixed_2020-09-22_123927.csv"
#equal <- "out/last_t_below_total_dyn_equal_2020-09-22_123617.csv"
#stable <- "out/last_t_below_total_dyn_stabilizing_2020-09-22_124640.csv"

file_list <- c(mixed, equal, stable)
make_plots <- function(competition){
  condition_type <- str_remove(competition, ".csv") %>% 
    str_remove("out/last_t_above_total_dyn_") %>% 
    str_split("_") %>% .[[1]] %>% .[1]
  
  last_t_out <- read_csv(competition)
  
  disp_rates <- sort(unique(last_t_out$dispersal))
  germ_fracs <- sort(unique(last_t_out$germination))
  surv_fracs <- sort(unique(last_t_out$survival))
  
  alpha_out <- last_t_out %>%
    group_by(dispersal, germination, survival, patch) %>% 
    mutate(N_pa = 1*(N > 0)) %>%
    summarize(alpha = sum(N_pa)) %>%
    ungroup()
  
  alpha_out <- alpha_out %>%
    group_by(dispersal, germination, survival) %>%
    summarize(mean_alpha = mean(alpha))
  
  gamma_out <- last_t_out %>%
    ungroup() %>%
    mutate(N_pa = 1*(N > 0)) %>%
    group_by(dispersal, germination, survival) %>%
    summarize(gamma = length(unique(species)))
  
  div_part <- alpha_out %>%
    left_join(gamma_out) %>%
    mutate(beta = gamma / mean_alpha)
  
  div_part %>%
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2))),
           survival = factor(survival, levels = surv_fracs, labels = paste("Surv. =", round(surv_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = partition, fill = partition)) +
    geom_point(alpha = 0.5) +
    geom_smooth() +
    facet_grid(rows = vars(germination), cols = vars(survival), drop = FALSE) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "top") +
    ggsave(paste0("figures/diversity_partitioning_", condition_type,".pdf"), width = 8, height = 12)
  
  div_part %>%
    mutate(#germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2))),
      survival = factor(survival, levels = surv_fracs, labels = paste("Surv. =", round(surv_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = germination, fill = germination, group = germination)) +
    geom_point(alpha = 0.25) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(partition), cols = vars(survival), drop = FALSE, scales = "free_y") +
    scale_x_log10() +
    #scale_y_log10() +
    scale_color_viridis_c("Germination", option = "B", end = .8) +
    scale_fill_viridis_c("Germination", option = "B", end = .8) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right") +
    ggsave(paste0("figures/germination_",condition_type,".pdf"), width = 12, height = 8)
  
  div_part %>%
    filter(germination <= 0.5) %>% 
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = survival, fill = survival, group = survival)) +
    geom_point(alpha = 0.25) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(partition), cols = vars(germination), drop = TRUE, scales = "free_y") +
    scale_x_log10() +
    #scale_y_log10() +
    scale_color_viridis_c("Survival", option = "D", end = .9) +
    scale_fill_viridis_c("Survival", option = "D", end = .9) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right") +
    ggsave(paste0("figures/survival_",condition_type,".pdf"), width = 12, height = 8)
}

### Loop through files
for(competition in file_list){
  make_plots(competition)
}


# anim<-div_part %>%
#   mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2)))) %>%
#   gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
#   mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"))) %>%
#   ggplot(aes(x = dispersal, y = diversity, color = survival, fill = survival, group = survival)) +
#   geom_point(alpha = 0.5) +
#   geom_smooth() +
#   facet_grid(rows = vars(partition), drop = FALSE, scales = "free_y") +
#   transition_manual(germination) +
#   enter_fade() +
#   exit_fade() +
#   scale_x_log10() +
#   #scale_y_log10() +
#   scale_color_viridis_c("Survival", option = "D", end = .9) +
#   scale_fill_viridis_c("Survival", option = "D", end = .9) +
#   labs(y = "Diversity", x = "Dispersal") +
#   theme(legend.position = "right")
# animate(anim,  renderer = ffmpeg_renderer(format = "mp4")) + 
#   anim_save("figures/test_vid.mp4")
