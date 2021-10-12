library("tidyverse")
library("broom")
library("colorspace")

normalize = function(x){
  normalized = (x-min(x))/(max(x)-min(x))
  return(normalized)
}

theme_set(theme_light() + 
            theme(strip.text = element_text(color = "black", size = 14),
                  strip.background = element_rect(fill = "gray90"),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  axis.ticks.length = unit(x = .2, units = "cm"),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 12)))

# read in raw data. Below, I'm using intermediate output, the code to create this is commented out but can be regenerated.
# sim_dat <- read_csv(file = "sim_output/final_sensitivity_2021-09-29_001433.csv")
# 
# # params <- sim_dat %>% 
# #   select(rep, dispersal, germination, survival, comp) %>% 
# #   unique()
# 
# nenvs <- 20
# nreps <- 250
# nspecies <- 40
# npatchs <- 100
# 
# nenvs * nreps * nspecies * npatchs * 2 == nrow(sim_dat)
# 
# # re-assign environment scenario to group by
# sim_dat$envscenario <- rep(1:nenvs, each = (nreps * nspecies * npatchs * 2))
# 
# alpha <- sim_dat %>% 
#   group_by(patch, envscenario, rep, dispersal, germination, survival, comp) %>% 
#   summarize(alpha = sum(1*(N>0))) %>%
#   ungroup() %>% 
#   group_by(rep, envscenario, comp) %>%
#   summarize(mean_alpha = mean(alpha))
# # left_join(params, by = c("rep", "comp"))
# write_csv(alpha, file = "summarized_output/sens_analysis_alphadiv.csv")
# 
# gamma <- sim_dat %>% 
#   mutate(N_pa = 1*(N > 0)) %>%
#   group_by(comp, envscenario, rep, dispersal, germination, survival, species) %>%
#   
#   summarize(presences = sum(N_pa)) %>% 
#   summarize(gamma = sum(presences > 0))
# write_csv(gamma, file = "summarized_output/sens_analysis_gammadiv.csv")

alpha <- read_csv("summarized_output/sens_analysis_alphadiv.csv")
gamma <- read_csv("summarized_output/sens_analysis_gammadiv.csv")

div_part <- alpha %>% 
  left_join(gamma) %>% 
  mutate(beta = gamma / mean_alpha) %>% 
  mutate(beta = ifelse(is.nan(beta), 0, beta))

# normalize and transform variables

div_part$mean_alpha <- normalize(div_part$mean_alpha)
div_part$beta <- normalize(div_part$beta)
div_part$gamma <- normalize(div_part$gamma)
div_part$dispersal <- normalize(log(div_part$dispersal))

# alpha_mod_eq <- lm(mean_alpha ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
# summary(alpha_mod_eq)
# 
# beta_mod_eq <- lm(beta ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
# summary(beta_mod_eq)
# 
# gamma_mod_eq <- lm(gamma ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "equal"))
# summary(gamma_mod_eq)
# 
# alpha_mod_stab <- lm(mean_alpha ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
# summary(alpha_mod_stab)
# 
# beta_mod_stab <- lm(beta ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
# summary(beta_mod_stab)
# 
# gamma_mod_stab <- lm(gamma ~ normalize(log(dispersal)) * germination * survival, data = subset(div_part, comp == "stable"))
# summary(gamma_mod_stab)




# try with quadratic terms
alpha_mod_eq <- lm(mean_alpha ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(alpha_mod_eq)

beta_mod_eq <- lm(beta ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(beta_mod_eq)

gamma_mod_eq <- lm(gamma ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(gamma_mod_eq)

alpha_mod_stab <- lm(mean_alpha ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
summary(alpha_mod_stab)

beta_mod_stab <- lm(beta ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
summary(beta_mod_stab)

gamma_mod_stab <- lm(gamma ~ dispersal + germination + survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
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

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("Intercept", "dispersal", "germination", "survival",
                            "I(dispersal^2)", "I(germination^2)", "I(survival^2)"))
sens_out <- sens_out %>% 
  mutate(term = recode_factor(term,
                              `Intercept` = "intercept",
                              `dispersal` = "dispersal",
                              `I(dispersal^2)` = "dispersal^2",
                              
                              `germination` = "germination",
                              `I(germination^2)` = "germination^2",
                              
                              `survival` = "survival",
                              `I(survival^2)` = "survival^2"
  ))

x_axis_labs <- c(expression(intercept), 
                 "dispersal",
                 expression(dispersal^2),
                 "germination",
                 expression(germination^2),
                 "survival",
                 expression(survival^2))

sensitivity_plot <- sens_out %>% 
  # mutate(term = factor(term, levels = c(
  #   "Intercept", "dispersal", "germination", "survival",
  #   "dispersal:germination", "dispersal:survival", "germination:survival",
  #   "dispersal:germination:survival"
  # ), labels = c(
  #   "intercept", "dispersal", "germination", "survival",
  #   "dispersal x \n germination", "dispersal x \n survival", "germination x \n survival",
  #   "dispersal x \n germination x \n survival"
  # ))) %>%
  
  ggplot(aes(y = estimate, ymin = (estimate - std.error), ymax = (estimate + std.error),
             color = response, fill = response, x = term)) +
  geom_bar(stat = "identity",
           position = position_dodge(), alpha = 0.7) +
  geom_vline(xintercept = c(1:6+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_errorbar(position = position_dodge(.9), width = 0.5, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.93,.62),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        axis.text.x = element_text(margin = margin(t = 10), vjust = 0)) +
  facet_wrap(~comp, nrow = 2) +
  labs(x = "", y = "Effect Size", fill = "Diversity level", color = "Diversity level") +
  scale_x_discrete(labels = x_axis_labs)

sensitivity_plot

ggsave(filename = "figures/sensitivity_analysis.pdf",plot = sensitivity_plot, width = 10, height = 8)

############################
# Incorporate interaction terms

alpha_mod_eq <- lm(mean_alpha ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "equal"))
summary(alpha_mod_eq)

beta_mod_eq <- lm(beta ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "equal"))
summary(beta_mod_eq)

gamma_mod_eq <- lm(gamma ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "equal"))
summary(gamma_mod_eq)

alpha_mod_stab <- lm(mean_alpha ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "stable"))
summary(alpha_mod_stab)

beta_mod_stab <- lm(beta ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "stable"))
summary(beta_mod_stab)

gamma_mod_stab <- lm(gamma ~ dispersal * germination * survival * I(dispersal^2) * I(germination^2) * I(survival^2), data = subset(div_part, comp == "stable"))
summary(gamma_mod_stab)

### 
alpha_mod_eq <- lm(mean_alpha ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(alpha_mod_eq)

beta_mod_eq <- lm(beta ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(beta_mod_eq)

gamma_mod_eq <- lm(gamma ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "equal"))
summary(gamma_mod_eq)

alpha_mod_stab <- lm(mean_alpha ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
summary(alpha_mod_stab)

beta_mod_stab <- lm(beta ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
summary(beta_mod_stab)

gamma_mod_stab <- lm(gamma ~ dispersal * germination * survival + I(dispersal^2) + I(germination^2) + I(survival^2), data = subset(div_part, comp == "stable"))
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

my_cols <- qualitative_hcl(n = 3)

sens_out$term <- factor(sens_out$term, 
                        levels = 
                          c("(Intercept)", "dispersal", "I(dispersal^2)", 
                            "germination", "I(germination^2)",
                            "survival", "I(survival^2)",
                            "dispersal:germination", "dispersal:survival", 
                            "germination:survival", "dispersal:germination:survival"))
x_axis_labs <- c("intercept", 
                 "dispersal",
                 expression(dispersal^2),
                 "germination",
                 expression(germination^2),
                 "survival",
                 expression(survival^2),
                 "dispersal:germination", "dispersal:survival", 
                 "germination:survival", "dispersal:germination:survival")

int_sensitivity_plot <- sens_out %>% 
  # mutate(term = factor(term, levels = c(
  #   "Intercept", "dispersal", "germination", "survival",
  #   "dispersal:germination", "dispersal:survival", "germination:survival",
  #   "dispersal:germination:survival"
  # ), labels = c(
  #   "intercept", "dispersal", "germination", "survival",
  #   "dispersal x \n germination", "dispersal x \n survival", "germination x \n survival",
  #   "dispersal x \n germination x \n survival"
  # ))) %>%
  
  ggplot(aes(y = estimate, ymin = (estimate - std.error), ymax = (estimate + std.error),
             color = response, fill = response, x = term)) +
  geom_vline(xintercept = c(1:10+0.5), alpha = 0.5, color = "gray50", size = 0.2) +
  geom_bar(stat = "identity",
           position = position_dodge(), alpha = 0.7) +
  geom_errorbar(position = position_dodge(.9), width = 0.5, show.legend = FALSE) +
  scale_color_manual(values = darken(my_cols, amount = .2)) + 
  scale_fill_manual(values = my_cols) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = c(.05,.92),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1)) +
  facet_wrap(~comp, nrow = 2) +
  labs(x = "", y = "Effect Size", fill = "Diversity level", color = "Diversity level") +
  scale_x_discrete(labels = x_axis_labs)

int_sensitivity_plot
ggsave("figures/sensitivity_analysis_int_quad.pdf", width = 15, height = 10)
