#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=ks_arms
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/ks_arms_%A_%a.log
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --array=0-11
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mck2199@columbia.edu

# One array task per (regional arm, anchor) of regional/ks_regional_arms.ipynb: 6 arms x 2 anchors.
# Each task executes its own copy of the notebook into logs/ and writes only its own caches.
# The combined tables come from one more run with FINAL=1 (all caches load, nothing recomputed):
#   jid=$(sbatch --parsable run/run_regional_arms_array.sbatch)
#   sbatch --array=0 --export=ALL,FINAL=1 --dependency=afterok:$jid run/run_regional_arms_array.sbatch
cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
export SSL_CERT_FILE=$(/burg-archive/home/mck2199/summer/.venv312/bin/python -c 'import certifi; print(certifi.where())')
export LD_LIBRARY_PATH=/burg-archive/home/mck2199/summer/.esmf/lib:$LD_LIBRARY_PATH
export ESMFMKFILE=/burg-archive/home/mck2199/summer/.esmf/lib/esmf.mk
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
# concurrent kernels sharing ~/.ipython corrupt its history sqlite; give each task its own
export IPYTHONDIR=$(mktemp -d)
trap 'rm -rf "$IPYTHONDIR"' EXIT

if [ -n "$FINAL" ]; then
  unset ARM ANCHOR
  NAME=all
else
  ARMS=(sas afr nae asia safca sasca)
  ANCHORS=(370 585)
  export ARM=${ARMS[$((SLURM_ARRAY_TASK_ID / 2))]}
  export ANCHOR=${ANCHORS[$((SLURM_ARRAY_TASK_ID % 2))]}
  NAME=${ARM}_${ANCHOR}
fi
echo "task $SLURM_ARRAY_TASK_ID: ${NAME}"

/burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert \
  --to notebook --execute \
  --ExecutePreprocessor.kernel_name=summer \
  --ExecutePreprocessor.timeout=21600 \
  --output-dir logs --output "ks_regional_arms_${NAME}.ipynb" \
  regional/ks_regional_arms.ipynb
