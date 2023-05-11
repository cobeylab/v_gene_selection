#!/bin/bash

simulations_dir=$1

# Concatenate files directory by directory, to avoid "argument list too long"
get_header=T

for dir in $simulations_dir/raw_simulation_files/*
do
    if [ $get_header == T ]
    then
        # For the first directory only, initialize combined file with header
        head $dir/simulation_individual* | head -2 | tail -1 > $simulations_dir/combined_simulations.csv
        get_header=F
    fi
    
    # Combine csv files with results
    tail -q -n+2 $dir/simulation_individual* >> $simulations_dir/combined_simulations.csv
    
    # Remove individual files
    #rm $dir/simulation_individual*

done

# Do the same for model parameter files
head $simulations_dir/raw_simulation_files/*/*model_parameters* | grep -v '==>' | head -1 > $simulations_dir/combined_model_parameters.csv
tail -q -n+2 $simulations_dir/raw_simulation_files/*/*model_parameters* >> $simulations_dir/combined_model_parameters.csv
