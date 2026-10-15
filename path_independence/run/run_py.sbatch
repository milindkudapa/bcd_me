#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=py
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/py_%j.log
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mck2199@columbia.edu

# Run one python script (argument) from path_independence/ with the ESMF environment, e.g.
#   sbatch --export=ALL,NULL_VARS=tasmax,tasmin,dtr run/run_py.sbatch method/audit/null1deg.py
cd /burg-archive/home/mck2199/summer/bcd_me/path_independence
source ~/.bashrc
export SSL_CERT_FILE=$(/burg-archive/home/mck2199/summer/.venv312/bin/python -c 'import certifi; print(certifi.where())')
export LD_LIBRARY_PATH=/burg-archive/home/mck2199/summer/.esmf/lib:$LD_LIBRARY_PATH
export ESMFMKFILE=/burg-archive/home/mck2199/summer/.esmf/lib/esmf.mk
export NUMBA_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "running $1 ($(date))"
/burg-archive/home/mck2199/summer/.venv312/bin/python "$1"
echo "exit $? ($(date))"
