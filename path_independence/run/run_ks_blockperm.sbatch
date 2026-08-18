#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=ks_blockperm
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/ks_%j.log
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

# Notebook to run: `sbatch run/run_ks_blockperm.sbatch <notebook.ipynb>`. A bare filename is
# looked up in the subfolders, so both `pipeline/ramip_summary.ipynb` and `ramip_summary.ipynb`
# work. cwd stays the package root because get_params() reads the dir_list.csv symlink there.
NB=${1:-pipeline/ramip_summary.ipynb}
if [ ! -f "$NB" ]; then
  for d in pipeline method earlier; do
    [ -f "$d/$NB" ] && NB="$d/$NB" && break
  done
fi
if [ ! -f "$NB" ]; then
  echo "notebook not found: $1" >&2
  exit 1
fi
echo "running $NB"

/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=28800 \
  "$NB"
