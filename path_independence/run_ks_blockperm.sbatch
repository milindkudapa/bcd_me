#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=ks_blockperm
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/ks_%j.log
#SBATCH --time=08:00:00
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

/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=28800 \
  ks_path_independence_ramip_multi.ipynb
