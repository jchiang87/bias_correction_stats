#!/bin/bash
#SBATCH --partition=rubin
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3:00:00

cd /sdf/home/j/jchiang/Bias_modeling_tests/baseline_tests
source ./setup.sh
#bps submit bps_cpDark_bias_median_per_row.yaml
bps restart --id u/jchiang/bias_dark_isr_13162_median_per_row_fp/20220521T000411Z
