#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=diag
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/diag_%j.log
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
/burg-archive/home/mck2199/summer/.venv312/bin/python /local/claude-560099/-burg-archive-home-mck2199-summer/3e3719b8-7718-46f2-9929-7c3e28ad3077/scratchpad/diag_miroc.py
