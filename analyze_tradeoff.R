library("tidyverse")
library("patchwork")

theme_set(theme_light() + 
            theme(strip.text = element_text(color = "black", size = 14),
                  strip.background = element_rect(fill = "gray90"),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  axis.ticks.length = unit(x = .2, units = "cm"),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 12)))

sim_dat <- read_csv("sim_output/final_tradeoff.csv")


alpha_out <- sim_dat %>% 
  group_by(tradeoff_strength, comp, patch, rep) %>% 
  mutate(N_pa = 1*(N > 0)) %>% 
  summarize(alpha = sum(N_pa)) %>% 
  ungroup() %>% 
  group_by(tradeoff_strength, comp, rep) %>% 
  summarize(mean_alpha = mean(alpha))

gamma_out <- sim_dat %>% 
  mutate(N_pa = 1*(N > 0)) %>% 
  filter(N_pa > 0) %>% 
  group_by(tradeoff_strength, comp, rep) %>% 
  summarize(gamma = length(unique(species)))

div_part <- alpha_out %>% 
  left_join(gamma_out) %>% 
  mutate(beta = gamma / mean_alpha) %>% 
  mutate(tradeoff_strength = factor(tradeoff_strength, levels = c("none", "weak", "strong")))

div_part <- div_part %>% 
  rename(alpha = mean_alpha)

equal_plot <- div_part %>% 
  pivot_longer(cols = c("alpha", "beta", "gamma"), names_to = "partition", values_to = "diversity") %>% 
  filter(comp == "equal") %>% 
  ggplot(aes(x = tradeoff_strength, y = diversity, color = partition)) + 
  geom_jitter(alpha = 0.5) +
  geom_boxplot(alpha = 0.2) +
  facet_grid(partition~comp, scales = "free") +
  theme(legend.position = "none")

stable_plot <- div_part %>% 
  pivot_longer(cols = c("alpha", "beta", "gamma"), names_to = "partition", values_to = "diversity") %>% 
  filter(comp == "stable") %>% 
  ggplot(aes(x = tradeoff_strength, y = diversity, color = partition)) + 
  geom_jitter(alpha = 0.5) +
  geom_boxplot(alpha = 0.2) +
  facet_grid(partition~comp, scales = "free") +
  theme(legend.position = "none")

equal_plot + stable_plot + 
  ggsave("figures/diversity_tradeoffs.pdf", width = 4*2, height = 3*3)
