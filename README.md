Examining diversity in metacommunities with dormant seed banks.

[![DOI](https://zenodo.org/badge/288286093.svg)](https://zenodo.org/badge/latestdoi/288286093)


This repository contains the code to reproduce the manuscript: 

**Wisnoski NI, and LG Shoemaker. "Seed banks alter metacommunity diversity: the interactive effects of competition, dispersal, and dormancy"**

To reproduce the analysis and figures from the manuscript, run the following:

1. `1_mc_sim.R` = this file includes the main simulation across parameter space in the paper. It generates output files that are used in figures 2 and 3.
2. `sim_output/extract_sim_output.sh` = commands included here parse the simulation output and split it into files for above and belowground data for equal and stable interactions. Paths and filenames will be user-specific and need to be specified for your system.
3.  `2_mc_sim_sensitivity.R` = this file includes the simulation model for the sensitivity analysis. It generates output files that are used to make figure 4. 
4.  `3_make_diversity_plots.R` = this file makes figures 2 and 3 from the output of the simulations reformated by the `extract_sim_output.sh` script.
5.  `4_plot_sensitivity_analysis.R` = this script creates figure 4 from the output of the sensitivity analysis simulation script.
6.  `s1_mc_sim_tradeoff.R` = simulation model for the scenario where species exhibit interspecific trait variation in dispersal, germination, and survival. 
7.  `s2_analyze_tradeoff.R` = this script makes the tradeoff related figures included in the supplement.
8.  `s3_visualize_example_tradeoffs.R` = this script demonstrates the differences in none, weak, and strong tradeoff scenarios explored in the supplement. 

