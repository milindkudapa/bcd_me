#!/bin/bash
#SBATCH --account=glab
#SBATCH --job-name=esmf_env
#SBATCH --output=/burg-archive/home/mck2199/summer/bcd_me/path_independence/logs/esmf_env_%j.log
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -e
CONDA=/burg/opt/anaconda3-2023.09/bin/conda
ENV=/burg-archive/home/mck2199/summer/.esmf
# nompi ESMF build: no MPI runtime needed at import time
[ -d $ENV ] || $CONDA create -y -p $ENV -c conda-forge --override-channels \
    python=3.12 "esmf=*=nompi*" esmpy

# esmpy is a pure-python wrapper; link it into the working venv and point it at the lib
SP=/burg-archive/home/mck2199/summer/.venv312/lib/python3.12/site-packages
ln -sfn $ENV/lib/python3.12/site-packages/esmpy $SP/esmpy
export LD_LIBRARY_PATH=/burg-archive/home/mck2199/summer/.esmf/lib:$LD_LIBRARY_PATH
export ESMFMKFILE=$ENV/lib/esmf.mk
/burg-archive/home/mck2199/summer/.venv312/bin/python - <<'PY'
import xesmf as xe, numpy as np, xarray as xr
ds = xr.Dataset(coords={'lat': np.arange(-89, 90, 2.0), 'lon': np.arange(0, 360, 2.0)})
out = xr.Dataset(coords={'lat': np.arange(-89.5, 90, 1.0), 'lon': np.arange(0.5, 360, 1.0)})
r = xe.Regridder(ds, out, 'bilinear', periodic=True)
print('xesmf OK', xe.__version__, r)
PY
