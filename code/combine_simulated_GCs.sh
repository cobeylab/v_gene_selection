#!/bin/bash

simulations_dir=$1

# Get csv header 
head -1 $simulations_dir/*simulation_individual* | head -2 | tail -1 > $simulations_dir/combined_simulations.csv

# Combine files
tail -q -n+2 $simulations_dir/simulation_individual_* >> $simulations_dir/combined_simulations.csv

# Remove originals

rm $simulations_dir/simulation_individual_*