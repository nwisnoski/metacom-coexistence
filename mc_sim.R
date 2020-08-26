library("tidyverse")
library("data.table")
library("som.nn")
library("progress")
library("RandomFields")


# define parameters

x_dim <- 100
y_dim <- 100
patches <- 100
species <- 20

# int mat params
intra = 1 
min_inter = 0 
max_inter = 1.5 
comp_scaler = 0.05

timesteps <- 1000
initialization <- 20
burn_in <- 80

# run sim

landscape <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)
env_df <- env_generate(landscape = landscape, env1Scale = 500, timesteps = timesteps+burn_in, plot = TRUE)

dynamics_out <- data.table()

# init_community
species_traits <- init_species(species)
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
  N_hat <- N*r / (1 + N%*%int_mat)
  N_hat[N_hat < 0] <- 0
  N_hat <- matrix(rpois(n = species*patches, lambda = N_hat), ncol = species, nrow = patches)
  
  # determine emigrants
  E <- matrix(rbinom(n = patches * species, size = N_hat, prob = rep(0.1, each = patches)))
  pb$tick()
}



