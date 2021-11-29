#!/bin/bash


# This file processes the simulation output to analyze in R

# first extract header row, then loop through sim runs for same session, 
# append individual files after removing first line with sed
head -1 total_dyn_equal_2021-03-15_143048.csv > final_equal.csv
for filename in $(ls total_dyn_equal_2021-03-1*.csv); do sed 1d $filename >> final_equal.csv; done
for filename in $(ls total_dyn_equal_2021-11-1*.csv); do sed 1d $filename >> final_equal.csv; done
head -1 total_dyn_stable_2021-03-15_222843.csv > final_stable.csv
for filename in $(ls total_dyn_stable_2021-03-1*.csv); do sed 1d $filename >> final_stable.csv; done
for filename in $(ls total_dyn_stable_2021-11-1*.csv); do sed 1d $filename >> final_stable.csv; done

# check to see that the header rows and columns look ok
head final_stable.csv
head final_equal.csv

# then split into above and belowground using awk to test for cols equal to 0
awk '$1!=0' FS=, final_equal.csv > final_equal_above.csv
awk '$2!=0' FS=, final_equal.csv > final_equal_below.csv
awk '$2!=0' FS=, final_stable.csv > final_stable_below.csv
awk '$1!=0' FS=, final_stable.csv > final_stable_above.csv
