#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=ks_1deg
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/ks1deg_%j.log
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mck2199@columbia.edu

cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
export SSL_CERT_FILE=$(/burg-archive/home/mck2199/summer/.venv312/bin/python -c 'import certifi; print(certifi.where())')
# ESMF from the conda side-env (built by build_esmf_env.sbatch)
export LD_LIBRARY_PATH=/burg-archive/home/mck2199/summer/.esmf/lib:$LD_LIBRARY_PATH
export ESMFMKFILE=/burg-archive/home/mck2199/summer/.esmf/lib/esmf.mk

/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=86400 \
  pipeline/ks_gwl_1deg.ipynb
