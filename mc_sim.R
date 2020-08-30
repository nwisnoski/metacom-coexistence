library("tidyverse")
library("data.table")
library("som.nn")
library("progress")
library("RandomFields")
source("metacom_functions.R")


# define parameters

x_dim <- 100
y_dim <- 100
patches <- 100
species <- 20
extirp_prob <- 0

# int mat params
intra = 1 
min_inter = 0 
max_inter = .5 
comp_scaler = 0.05

# seed bank params
#surv <- 0.8 # survival rate in the seed bank
#germ <- 0.5 # fraction of seeds that germinate
#zeta <- 0 # magnitude of the effect of env. stochasticity


timesteps <- 2200
initialization <- 200
burn_in <- 800

# run sim

landscape <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)
env_df <- env_generate(landscape = landscape, env1Scale = 500, timesteps = timesteps+burn_in, plot = TRUE)

dynamics_out <- data.table()

# init_community
species_traits <- init_species(species, dispersal_rate = 0.001, env_niche_breadth = 0.25)
disp_array <- generate_dispersal_matrices(landscape, species, patches, species_traits, torus = FALSE)
int_mat <- species_int_mat(species = species, intra = intra, 
                           min_inter = min_inter, max_inter = max_inter,
                           comp_scaler = comp_scaler, plot = TRUE)

pb <- progress_bar$new(total = initialization + burn_in + timesteps)
N <- init_community(initialization = initialization, species = species, patches = patches)

for(i in 1:(initialization + burn_in + timesteps)){
  if(i <= initialization){
    if(i %in% seq(10, 100, by = 10)){
      N <- N + matrix(rpois(n = species*patches, lambda = 0.5), nrow = patches, ncol = species)
    }
    env <- env_df$env1[env_df$time == 1]
  } else {
    env <- env_df$env1[env_df$time == (i - initialization)]
  }
  
  # compute r
  r <- compute_r_xt(species_traits, env = env, species = species)
  
  # compute growth
  #N_hat <- N*r / (1 + N%*%int_mat)
  N_hat <- survival(N, species_traits) + growth(N, species_traits, r, int_mat)
  N_hat[N_hat < 0] <- 0
  N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches)
  
  # determine emigrants
  E <- matrix(nrow = patches, ncol = species)
  for(s in 1:species){
    E[,s] <- rbinom(n = patches, size = N_hat[,s], prob = species_traits$dispersal_rate[s])
  }
  dispSP <- colSums(E)
  
  I_hat_raw <- disp_array[,,1]%*%E
  I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
  I_hat[is.nan(I_hat)] <- 1
  I <- sapply(1:species, function(x) {
    if(dispSP[x]>0){
      table(factor(sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]), levels = 1:patches))
    } else {rep(0, patches)}
  })
  
  
  I_hat_raw <- matrix(nrow = patches, ncol = species)
  
  for(s in 1:species){
    I_hat_raw[,s] <- disp_array[,,s] %*% E[,s]
  }
  
  # standardize so colsums = 1
  I_hat <- t(t(I_hat_raw)/colSums(I_hat_raw))
  I_hat[is.nan(I_hat)] <- 1
  
  I <- sapply(1:species, function(x) {
    if(dispSP[x]>0){
      table(factor(
        sample(x = patches, size = dispSP[x], replace = TRUE, prob = I_hat[,x]),
        levels = 1:patches))
    } else {rep(0, patches)}
  })
  
  N <- N_hat - E + I
  
  N[rbinom(n = species * patches, size = 1, prob = extirp_prob) > 0] <- 0
  
  dynamics_out <- rbind(dynamics_out, 
                        data.table(N = c(N),
                                   patch = 1:patches,
                                   species = rep(1:species, each = patches),
                                   env = env,
                                   time = i-initialization-burn_in))
  
  pb$tick()
}

dynamics_out %>% 
  filter(time %in% seq(min(dynamics_out$time), max(dynamics_out$time), by = 10)) %>% 
  ggplot(aes(x = time, y = N, group = species, color = species)) + 
  geom_line(alpha = 0.3) + 
  facet_wrap(~patch) +
  scale_color_viridis_c() +
  theme_minimal()

end_sbs <- dynamics_out %>% 
  filter(time == max(dynamics_out$time)) %>% 
  select(-env, -time) %>% 
  spread(key = species, value = N) %>% 
  select(-patch)

library(vegan)
vegan::
adipart(end_sbs)

