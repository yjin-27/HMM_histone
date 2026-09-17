#!/bin/bash
#SBATCH -A b1042
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --partition  genomics
#SBATCH --job-name=gaussian_pyhmm_2-7states_H3K4me3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuchengjin2027@u.northwestern.edu

# Adapted from Liebe's pyhhmm.sh.
# Trains 2-7 state Gaussian HMMs on raw coverage of the two ENCSR000DWD
# isogenic replicates (chr14 by default; export CHROM=chrN to override).
# Outputs per replicate x n_states:
#   *parameters.txt   (pi, emission means/variances, transition matrix, logL history)
#   *regions.bed      (colored state segmentation -- load this in IGV)
#   *transitions.bed  (state transition points)

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
"$PYBIN" gaussian_hmm_histone.py
