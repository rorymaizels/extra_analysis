#!/bin/bash
#SBATCH --job-name=length
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=250G
#SBATCH --partition=cpu
#SBATCH --output=len-%j.out

script=length_analysis.py
environment=velvet_env1

python_path=/camp/home/maizelr/.conda2/my_envs/$environment/bin
$python_path/python $script