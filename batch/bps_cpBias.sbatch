#!/bin/bash
#SBATCH --partition=rubin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=30:00

cd /sdf/home/j/jchiang/Bias_modeling_tests/baseline_tests
source ./setup.sh
#bps submit bps_cpBias_median_per_row.yaml
bps restart --id u/jchiang/bias_13162_median_per_row/20220520T213227Z
