library("tidyverse")
library("data.table")
library("som.nn")
library("progress")
library("RandomFields")
library("vegan")
library("foreach")
library("doParallel")
source("metacom_functions.R")

dynamics_total <- read_csv("sim_output/total_dyn_mixed_2020-09-16_150314.csv")
# analyze diversity
last_t_out <- dynamics_total %>% 
  filter(time == max(dynamics_total$time),
         N > 0) %>% 
  group_by(dispersal, germination, survival, patch) %>% 
  arrange(patch)

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
  theme(legend.position = "top") +
  ggsave("figures/aboveground_diversity_partitioning_mixed2000.pdf", width = 8, height = 6)


last_t_out <- dynamics_total %>% 
  filter(time == max(dynamics_total$time),
         D > 0) %>% 
  group_by(dispersal, germination, survival, patch) %>% 
  arrange(patch)

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
  theme(legend.position = "top") +
  ggsave("figures/belowground_diversity_partitioning_mixed2000.pdf", width = 8, height = 6)

