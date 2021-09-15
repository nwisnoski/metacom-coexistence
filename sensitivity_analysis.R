library("tidyverse")
library("broom")
library("colorspace")

normalize = function(x){
  normalized = (x-min(x))/(max(x)-min(x))
  return(normalized)
}


sim_dat <- read_csv(file = "sim_output/final_sensitivity_2021-09-13_123137.csv")

params <- sim_dat %>% 
  select(rep, dispersal, germination, survival, comp) %>% 
  unique()

alpha <- sim_dat %>% 
  group_by(patch, rep, dispersal, germination, survival, comp) %>% 
  mutate(N_pa = 1*(N > 0)) %>%
  summarize(alpha = sum(N_pa)) %>%
  ungroup() %>% 
  group_by(rep, comp) %>%
  summarize(mean_alpha = mean(alpha)) %>% 
  left_join(params, by = c("rep", "comp"))

gamma <- sim_dat %>% 
  mutate(N_pa = 1*(N > 0)) %>%
  group_by(comp, rep, dispersal, germination, survival, species) %>%

  summarize(presences = sum(N_pa)) %>% 
  summarize(gamma = sum(presences > 0))

div_part <- alpha %>% 
  left_join(gamma) %>% 
  mutate(beta = gamma / mean_alpha)

alpha_mod_eq <- lm(mean_alpha ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
summary(alpha_mod_eq)

beta_mod_eq <- lm(beta ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
summary(beta_mod_eq)

gamma_mod_eq <- lm(gamma ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
summary(gamma_mod_eq)

alpha_mod_stab <- lm(mean_alpha ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
summary(alpha_mod_stab)

beta_mod_stab <- lm(beta ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
summary(beta_mod_stab)

gamma_mod_stab <- lm(gamma ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
summary(gamma_mod_stab)


# build dataframe to plot
sens_out <- tidy(alpha_mod_eq) %>% 
  mutate(response = "alpha", comp = "equal") %>% 
  bind_rows(
    (tidy(beta_mod_eq) %>% 
       mutate(response = "beta", comp = "equal")), 
    (tidy(gamma_mod_eq) %>% 
       mutate(response = "gamma", comp = "equal")),
    
    (tidy(alpha_mod_stab) %>% 
       mutate(response = "alpha", comp = "stable")),
    (tidy(beta_mod_stab) %>% 
       mutate(response = "beta", comp = "stable")), 
    (tidy(gamma_mod_stab) %>% 
       mutate(response = "gamma", comp = "stable"))
      ) 

sens_out$term <- str_replace(string = sens_out$term, pattern = coll("(Intercept)"), replacement =  "Intercept")
sens_out$term <- str_wrap(str_replace(string = sens_out$term, pattern = coll("normalize(log(dispersal))"), replacement =  "dispersal"), width = 10)


my_cols <- qualitative_hcl(n = 3)

sensitivity_plot <- sens_out %>% 
  mutate(term = factor(term, levels = c(
    "Intercept", "dispersal", "germination", "survival",
    "dispersal:germination", "dispersal:survival", "germination:survival",
    "dispersal:germination:survival"
  ), labels = c(
    "intercept", "dispersal", "germination", "survival",
    "dispersal x \n germination", "dispersal x \n survival", "germination x \n survival",
    "dispersal x \n germination x \n survival"
  ))) %>% 
  ggplot(aes(y = estimate, ymin = (estimate - std.error), ymax = (estimate + std.error),
             color = response, fill = response, x = term)) +
  geom_bar(stat = "identity",
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(.9), width = 0.5, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  theme_bw() + theme(panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     legend.position = c(.08,.9)) +
  facet_wrap(~comp, nrow = 2) +
  labs(x = "", y = "Effect Size", fill = "Diversity level", color = "Diversity level")

sensitivity_plot

ggsave(filename = "figures/sensitivity_analysis.pdf",plot = sensitivity_plot, width = 10, height = 8)
  
  
