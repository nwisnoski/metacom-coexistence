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

sim_dat <- read_csv("sim_output/final_tradeoff_2021-07-20_172050.csv")


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
  rename(alpha = mean_alpha) %>%
  pivot_longer(cols = c("alpha", "beta", "gamma"), names_to = "partition", values_to = "diversity")

div_part_meansd <- div_part %>% 
  group_by(tradeoff_strength, comp, partition) %>% 
  summarize(mean_diversity = mean(diversity),
            sd_diversity = sd(diversity)) %>% 
  right_join(div_part)

equal_plot <- div_part_meansd %>% 
  #pivot_longer(cols = c("alpha", "beta", "gamma"), names_to = "partition", values_to = "diversity") %>% 
  filter(comp == "equal") %>% 
  ggplot(aes(x = tradeoff_strength, y = diversity)) + 
  geom_jitter(alpha = 0.5, width = 0.2) +
  geom_point(aes(y = mean_diversity), size = 3) +
  geom_errorbar(aes(ymax = mean_diversity + sd_diversity, 
                    ymin = mean_diversity - sd_diversity), 
                size = 1, 
                width = 0.2) + 
  facet_grid(partition~comp, scales = "free") +
  theme(legend.position = "none") + 
  labs(x = "Trade-off strength", y = "Diversity")

stable_plot <- div_part_meansd %>% 
  #pivot_longer(cols = c("alpha", "beta", "gamma"), names_to = "partition", values_to = "diversity") %>% 
  filter(comp == "stable") %>% 
  ggplot(aes(x = tradeoff_strength, y = diversity)) + 
  geom_jitter(alpha = 0.5, width = 0.2) +
  geom_point(aes(y = mean_diversity), size = 3) +
  geom_errorbar(aes(ymax = mean_diversity + sd_diversity, 
                    ymin = mean_diversity - sd_diversity), 
                size = 1, 
                width = 0.2) + 
  facet_grid(partition~comp, scales = "free") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(x = "Trade-off strength", y = "Diversity")

equal_plot + stable_plot + 
  ggsave("figures/diversity_tradeoffs.pdf", width = 4.5*2, height = 3*3)


# Scales on y-axis are hard to compare. So that's why I made two figures above and joined
div_part_meansd %>% 
  ggplot(aes(x = tradeoff_strength, y = diversity)) + 
  geom_jitter(alpha = 0.5, width = 0.2, height = 0) +
  geom_point(aes(y = mean_diversity), size = 3) +
  geom_errorbar(aes(ymax = mean_diversity + sd_diversity, 
                    ymin = mean_diversity - sd_diversity), 
                size = 1, 
                width = 0.2) + 
  facet_grid(partition~comp, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Trade-off strength", y = "Diversity")

#####
sim_dat %>% 
  group_by(rep) %>% 
  ggplot(aes(x = dispersal, y = germination, color = rep)) + 
  geom_point() +
  facet_grid(comp ~ tradeoff_strength) +
  theme(axis.text = element_text(size = 10))



### Quant analysis of traits
trait_correlations <- sim_dat %>% 
  group_by(tradeoff_strength, comp, rep) %>% 
  summarize(disp_germ_corr = cor(dispersal, germination),
            disp_surv_corr = cor(dispersal, survival),
            germ_surv_corr = cor(germination, survival))

ddisp_germ_equal_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "equal") %>% 
  ggplot(aes(x = disp_germ_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none") +
  labs(title = "Dispersal-Germination", x = "Trait correlation", y = "Diversity")
disp_germ_stable_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "stable") %>% 
  ggplot(aes(x = disp_germ_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none") +
  labs(title = "Dispersal-Germination", x = "Trait correlation", y = "Diversity")
disp_surv_equal_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "equal") %>% 
  ggplot(aes(x = disp_surv_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none", 
        axis.title.y = element_blank()) +
  labs(title = "Dispersal-Survival", x = "Trait correlation", y = "Diversity")
disp_surv_stable_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "stable") %>% 
  ggplot(aes(x = disp_surv_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(title = "Dispersal-Survival", x = "Trait correlation", y = "Diversity")
germ_surv_equal_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "equal") %>% 
  ggplot(aes(x = germ_surv_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(title = "Germination-Survival", x = "Trait correlation", y = "Diversity")
germ_surv_stable_plot <- div_part_meansd %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep")) %>% 
  filter(comp == "stable") %>% 
  ggplot(aes(x = germ_surv_corr, y = diversity)) + 
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm") +
  facet_grid(partition~., scales = "free") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  labs(title = "Germination-Survival", x = "Trait correlation", y = "Diversity")





disp_germ_stable_plot + disp_surv_stable_plot + germ_surv_stable_plot +
  ggsave("figures/diversity_tradeoffs_quantitative_stable.pdf",  width = 5*2, height = 2*3)
disp_germ_equal_plot + disp_surv_equal_plot + germ_surv_equal_plot +
  ggsave("figures/diversity_tradeoffs_quantitative_equal.pdf",  width = 5*2, height = 2*3)
