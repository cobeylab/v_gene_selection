#!/bin/bash

simulations_dir=$1

# Get csv header 
head $simulations_dir/raw_simulation_files/*/simulation_individual* | head -2 | tail -1 > $simulations_dir/combined_simulations.csv

# Combine csv files with results
tail -q -n+2 $simulations_dir/raw_simulation_files/*/simulation_individual* >> $simulations_dir/combined_simulations.csv

# Remove originals

rm $simulations_dir/raw_simulation_files/*/simulation_individual*

# Do the same for model parater files
head $simulations_dir/raw_simulation_files/*/*model_parameters* | head -2 | tail -1 > $simulations_dir/combined_model_parameters.csv
tail -q -n+2 $simulations_dir/raw_simulation_files/*/*model_parameters* >> $simulations_dir/combined_model_parameters.csv
