# metacom-coexistence
Examining coexistence in metacommunities with dormancy.

This repository contains the code to reproduce the manuscript: 

Wisnoski NI, and LG Shoemaker. "Seed banks alter metacommunity diversity: the interactive effects of competition, germination, and survival"

To reproduce the analysis and figures from the manuscript, run the following:

1. `mc_sim.R` = this file includes the main simulation across parameter space in the paper.
2. `sim_output/extract_sim_output.sh` = commands included here parse the simulation output and split it into files for above and belowground data for equal and stable interactions. Paths and filenames here may be user-specific and likely need tweaking for your system.
3.  `make_plots.R` = this file makes figures from the output of the simulations
4.  `mc_sim_tradeoff.R` = simulation for condition where there's interspecific trait variation
5.  `analyze_tradeoff.R` = making tradeoff related figures.

