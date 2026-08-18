#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=ks_empnull
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/ks_empnull_%j.log
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mck2199@columbia.edu

cd /burg-archive/home/mck2199/summer/bcd_me/path_independence

# ARRAYLAKE_TOKEN must reach the job (it lives in ~/.bashrc)
source ~/.bashrc

# regionmask/cartopy download Natural Earth shapefiles; cluster SSL needs certifi
export SSL_CERT_FILE=$(/burg-archive/home/mck2199/summer/.venv312/bin/python -c 'import certifi; print(certifi.where())')

export KS_NWORKERS=${SLURM_CPUS_PER_TASK:-16}

# 1) empirical-null analysis (new)
/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=12600 \
  earlier/ks_empirical_null.ipynb

# 2) re-render the pairwise figure with authors' cmap_mod (loads cache, no recompute)
/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=3600 \
  earlier/ks_path_independence_maps.ipynb
