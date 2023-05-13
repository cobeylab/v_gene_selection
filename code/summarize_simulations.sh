#!/bin/bash
results_dir=$1

sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=summarize_simulations
#SBATCH --output=out_err_files/summarize_simulations.out
#SBATCH --error=out_err_files/summarize_simulations.err
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=100:00:00
#SBATCH --mem-per-cpu=256G

module load R/3.6.1
Rscript summarize_simulations.R $results_dir

EOT
