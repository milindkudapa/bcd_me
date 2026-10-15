#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=dl_ssp585
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/dl_ssp585_%j.log
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mck2199@columbia.edu

# ssp585 tasmax + tasmin from pangeo, same members and layout as the tas/pr download,
# so the regional-arms notebooks can test them against the RAMIP arms.
cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
export SSL_CERT_FILE=$(/burg-archive/home/mck2199/summer/.venv312/bin/python -c 'import certifi; print(certifi.where())')
export IPYTHONDIR=$(mktemp -d)
trap 'rm -rf "$IPYTHONDIR"' EXIT
export DL_EXPS=ssp585 DL_VARS=tasmax,tasmin

/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=21600 \
  --output-dir logs --output "scenariomip_download_ssp585_tmax.ipynb" \
  pipeline/scenariomip_download.ipynb
