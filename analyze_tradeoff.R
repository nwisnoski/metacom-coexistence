library("tidyverse")
library("patchwork")
library("ggridges")

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





div_corrs <- div_part_meansd %>% 
  # filter(partition == "alpha",
  #        comp == "equal") %>% 
  left_join(trait_correlations, by = c("tradeoff_strength", "comp", "rep"))


trait_summaries <- sim_dat %>% 
  group_by(tradeoff_strength, comp) %>% 
  summarize(disp_mean = mean(dispersal),
            disp_sd = sd(dispersal),
            germ_mean = mean(germination),
            germ_sd = sd(germination),
            surv_mean = mean(survival),
            surv_sd = sd(survival)) %>% 
  group_by(tradeoff_strength, comp)

# This data wrangling tip came from https://community.rstudio.com/t/pivot-longer-on-multiple-column-sets-pairs/43958/10
trait_plot <- trait_summaries %>% 
  pivot_longer(cols = disp_mean:surv_sd, 
               names_to = c("trait", ".value"),
               names_pattern = "(.+)_(.+)") %>% 
  mutate(trait = case_when(
    trait == "disp" ~ "Dispersal",
    trait == "germ" ~ "Germination",
    trait == "surv" ~ "Survival"
  )) %>% 
  
  ggplot(aes(x = tradeoff_strength, y = mean, ymin = mean-sd, ymax = mean+sd, color = comp)) + 
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = .5)) + 
  geom_errorbar(width = .5, size = .5, position = position_dodge(width = .5)) +
  facet_grid(~trait, scales = "free") +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10)) +
  labs(y = "Trait mean (±s.d.)", x = "Correlation strength", 
       color = "Competition") +
  ggsave("figures/trait_values.pdf", width = 4*3, height = 3)


## Make figure 4
div_corr_plot <- div_part_meansd %>% 
  ggplot(aes(x = diversity, y = tradeoff_strength)) + 
  geom_density_ridges() + 
  facet_grid(partition~comp, scales = "free_x") +
  labs(x = "Diversity", y = "Trait correlation strength") +
  ggsave("figures/diversity_covariation.pdf", width = 8, height = 8*3/4)
