#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=gdp
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/gdp_%A_%a.log
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mck2199@columbia.edu

# Stage 1, the country-temperature cache, one array task per model:
#   sbatch --array=0-8 gdp/run_gdp.sbatch
# Stage 2, Burke's regressions (self-check) and the two notebooks, after the array:
#   sbatch --dependency=afterok:<array job id> --export=ALL,STAGE=nb gdp/run_gdp.sbatch
#   (NBS="gdp/ramip_gdp_pathway.ipynb" runs one notebook; several are separated by spaces)
cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
PY=/burg-archive/home/mck2199/summer/.venv312/bin/python
export SSL_CERT_FILE=$($PY -c 'import certifi; print(certifi.where())')
export LD_LIBRARY_PATH=/burg-archive/home/mck2199/summer/.esmf/lib:$LD_LIBRARY_PATH
export ESMFMKFILE=/burg-archive/home/mck2199/summer/.esmf/lib/esmf.mk
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
MODELS=(CanESM5-1 CESM2 CNRM-ESM2-1 EC-Earth3-AerChem GISS-E2-1-G MIROC6 MRI-ESM2-0 NorESM2-LM UKESM1-0-LL)

if [ "${STAGE:-cache}" = cache ]; then
  echo "cache ${MODELS[$SLURM_ARRAY_TASK_ID]} ($(date))"
  $PY gdp/country_tas.py "${MODELS[$SLURM_ARRAY_TASK_ID]}"
else
  $PY gdp/burke.py || exit 1
  export IPYTHONDIR=$(mktemp -d)                 # a private IPython history db per job
  for NB in ${NBS:-gdp/ramip_gdp_burke.ipynb gdp/ramip_gdp_pathway.ipynb}; do
    echo "running $NB ($(date))"
    /burg-archive/home/mck2199/summer/.venv312/bin/jupyter nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.kernel_name=summer --ExecutePreprocessor.timeout=10000 "$NB" || exit 1
  done
fi
echo "exit $? ($(date))"
