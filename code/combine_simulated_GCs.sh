#!/bin/bash

simulations_dir=$1

# Get csv header 
head $simulations_dir/raw_simulation_files/*/simulation_individual* | head -2 | tail -1 > $simulations_dir/combined_simulations.csv

# Combine files
tail -q -n+2 $simulations_dir/raw_simulation_files/*/simulation_individual* >> $simulations_dir/combined_simulations.csv

# Remove originals

rm $simulations_dir/simulation_individual_*