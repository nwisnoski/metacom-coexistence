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
equal_above <- "sim_output/final_equal_above.csv"
stable_above <- "sim_output/final_stable_above.csv"
equal_below <- "sim_output/final_equal_below.csv"
stable_below <- "sim_output/final_stable_below.csv"


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
    filter(survival %in% c(0.5, 0.99))
  
  
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
        survival == 0.5 ~ "Survival = 0.5",
        survival == 0.99 ~ "Survival = 0.99"
      )
    ) %>% bind_rows(no_sb)
  
  
  germ_plot <- div_part %>%
  
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"),
                              labels = c("mean alpha", "beta", "gamma"))) %>%
    ggplot(aes(x = dispersal, y = diversity, color = germination, fill = germination, group = germination)) +
    geom_point(alpha = 0.05) +
    #geom_line(mapping = aes(group = paste(germination, rep)), stat = "smooth", alpha = 0.2, se = F) +
    geom_smooth(alpha = 0.25) +
    facet_grid(rows = vars(partition), cols = vars(survival), drop = FALSE, scales = "free_y") +
    scale_x_log10() +
    scale_y_continuous(breaks = scales::pretty_breaks(4)) +
    scale_color_viridis_c("Germination", option = "B", end = .8) +
    scale_fill_viridis_c("Germination", option = "B", end = .8) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right")
  ggsave(filename = paste0("figures/germination_", position, "_", condition_type,".pdf"),
         plot = germ_plot, width = 12, height = 8)
  
  surv_plot <- div_part %>%
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
    scale_y_continuous(breaks = scales::pretty_breaks(4)) +
    scale_color_viridis_d("Survival", option = "D", end = .9) +
    scale_fill_viridis_d("Survival", option = "D", end = .9) +
    labs(y = "Diversity", x = "Dispersal") +
    theme(legend.position = "right") 
  ggsave(filename = paste0("figures/survival_", position, "_", condition_type,".pdf"), 
         plot = surv_plot, width = 12, height = 8)
  
  germ_comp <- div_part %>% 
    filter(germination == 1 | germination == 0.1) %>% 
    filter(survival == 1) %>% 
    mutate(germination = factor(germination, levels = germ_fracs, labels = paste("Germ. =", round(germ_fracs,2))),
           survival = factor(survival, levels = surv_fracs, labels = paste("Surv. =", round(surv_fracs,2)))) %>%
    gather(mean_alpha, beta, gamma, key = "partition", value = "diversity") %>%
    mutate(partition = factor(partition, levels = c("mean_alpha", "beta", "gamma"), 
                              labels = c("mean alpha", "beta", "gamma"))) 
  
  
}

### Loop through files
for(competition in file_list){
  make_plots(competition)
}



