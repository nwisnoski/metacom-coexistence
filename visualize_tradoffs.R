library("tidyverse")
library("patchwork")

species <- 40

tradeoff_strength <- c("strong", "weak", "none")


tradeoff_df <- data.frame()

for(y in tradeoff_strength){
  if(y == "strong"){
    tradeoff_noise = .25
    disp_tradeoff <- -0.5
    germ_tradeoff <- -0.5
    surv_tradeoff <- 0.5
  }
  if(y == "weak"){
    tradeoff_noise = 1
    disp_tradeoff <- -0.5
    germ_tradeoff <- -0.5
    surv_tradeoff <- 0.5
  }
  if(y == "none"){
    tradeoff_noise = 1
    disp_tradeoff = 0
    germ_tradeoff = 0
    surv_tradeoff = 0
  }
  
  # seed_mass = exp(rnorm(n = species, mean = 0, sd = 3.45/2))
  # sd(seed_mass)/mean(seed_mass)
  # hist(seed_mass,breaks = 30)
  # retrying with a normal distribution, too much variance in lognormal
  seed_mass <- rnorm(n = species, mean = 4.14, sd = 2)
  hist(seed_mass, breaks = 30)
  seed_mass <- seed_mass + abs(min(seed_mass))
  hist(seed_mass, breaks = 30)
  
  
  # then, based on known trade-offs, we generate the corresponding trait values
  disp <- seed_mass * disp_tradeoff + rnorm(species, sd = tradeoff_noise)
  disp <- (disp - min(disp)) / (max(disp) - min(disp)) *1/1000 # multiply by 1/1000 to range from 0 to 0.001
  plot(seed_mass, disp)
  
  germ <- seed_mass * germ_tradeoff + rnorm(species, sd = tradeoff_noise)
  germ <- (germ - min(germ)) / (max(germ) - min(germ))
  plot(seed_mass, germ)
  
  surv <- seed_mass * surv_tradeoff + rnorm(species, sd = tradeoff_noise)
  surv <- (surv - min(surv)) / (max(surv) - min(surv))
  plot(seed_mass, surv)
  
  tradeoffs <- data.frame(
    tradeoff = y,
    seed_mass = seed_mass,
    dispersal = disp, 
    germination = germ,
    survival = surv
  )
  
  tradeoff_df <- rbind.data.frame(tradeoff_df, tradeoffs)
  
}

germ_disp_fig <- tradeoff_df %>% 
  ggplot(aes(x = dispersal, y = germination)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tradeoff) +
  theme(axis.text = element_text(size=8)) +
  scale_x_continuous(labels = scales::scientific_format())

disp_surv_fig <- tradeoff_df %>% 
  ggplot(aes(x = dispersal, y = survival)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tradeoff) +
  theme(axis.text = element_text(size=8)) +
  scale_x_continuous(labels = scales::scientific_format())


germ_surv_fig <- tradeoff_df %>% 
  ggplot(aes(x = survival, y = germination)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~tradeoff) +
  theme(axis.text = element_text(size=8))


germ_disp_fig + disp_surv_fig + germ_surv_fig + 
  plot_layout(nrow = 3) +
  ggsave("figures/tradeoff_examples.pdf", width = 9, height = 9)
