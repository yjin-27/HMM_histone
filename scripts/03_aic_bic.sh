#!/bin/bash
#SBATCH -A b1042
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --partition  genomics
#SBATCH --job-name=aic_bic_H3K4me3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuchengjin2027@u.northwestern.edu

# Computes AIC/BIC across the 2-7 state models (replicating the
# "AIC BIC num states" section of Liebe's pyhhmm.ipynb) and writes a
# CSV table + bar plots. Light enough for a login node too:
#   conda activate hmmchip && python scripts/src/aic_bic_histone.py

module load python-miniconda3
source activate hmmchip

set -euo pipefail
unset PYTHONPATH PYTHONHOME

# Quest gotcha: 'module load python-miniconda3' prepends the system Python
# to PATH and can shadow the conda env even when (hmmchip) shows in the
# prompt. Resolve the env's own interpreter explicitly and use it directly.
PYBIN="$(conda env list | awk '$1=="hmmchip"{print $NF}')/bin/python"
[ -x "$PYBIN" ] || { echo "ERROR: hmmchip python not found at $PYBIN"; exit 1; }
echo "Using interpreter: $PYBIN"
"$PYBIN" -c "import pyBigWig" || { echo "ERROR: wrong interpreter (pyBigWig missing)"; exit 1; }

cd /projects/b1042/AmaralLab/yucheng/HMM_histone/scripts/src
"$PYBIN" aic_bic_histone.py
