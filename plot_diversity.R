# analyze diversity
library("tidyverse")
library("data.table")

plot_diversity <- function(dynamics_total, file){
  
  # plot aboveground
  last_t_out <- dynamics_total %>%
    filter(time == max(dynamics_total$time),
           N > 0) %>%
    group_by(dispersal, germination, survival, patch) %>%
    arrange(patch)
  write_csv(last_t_out, path = paste0("out/last_t_above_", file))
  
  alpha_out <- last_t_out %>%
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
    theme_light() +
    scale_x_log10() +
    scale_y_log10() + 
    labs(y = "Diversity") +
    theme(legend.position = "top") +
    ggsave(paste0("figures/diversity_partitioning_",str_remove(file,".csv"),"_above.pdf"), width = 8, height = 6)
  
  # plot seedbank
  last_t_out <- dynamics_total %>%
    filter(time == max(dynamics_total$time),
           D > 0) %>%
    group_by(dispersal, germination, survival, patch) %>%
    arrange(patch)
  write_csv(last_t_out, path = paste0("out/last_t_below_", file))
  
  alpha_out <- last_t_out %>%
    mutate(D_pa = 1*(D > 0)) %>%
    summarize(alpha = sum(D_pa)) %>%
    ungroup()
  
  alpha_out <- alpha_out %>%
    group_by(dispersal, germination, survival) %>%
    summarize(mean_alpha = mean(alpha))
  
  gamma_out <- last_t_out %>%
    ungroup() %>%
    mutate(D_pa = 1*(D > 0)) %>%
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
    theme_light() +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Diversity") +
    theme(legend.position = "top") +
    ggsave(paste0("figures/diversity_partitioning_",str_remove(file,".csv"),"_below.pdf"), width = 8, height = 6)
  
}

disp_rates <- 10^seq(-5, 0, length.out = 50)
germ_fracs <- seq(.1,1, length.out = 4)
surv_fracs <- seq(.1,1, length.out = 4)


files <- list.files("sim_output/stochastic/")
for(f in files){
  if(!(file.exists(paste0("figures/diversity_partitioning_",str_remove(f,".csv"),"_above.pdf")) & 
       file.exists(paste0("figures/diversity_partitioning_",str_remove(f,".csv"),"_below.pdf")))){
    dat <- data.table::fread(paste0("sim_output/stochastic/",f))
    plot_diversity(dat, file = f)
  }
}
