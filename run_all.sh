#!/bin/bash
# ============================================================
# ONE-COMMAND DRIVER for the H3K4me3 HMM workflow.
#
# Usage (from a Quest login node, ideally inside tmux):
#
#     cd /projects/b1042/AmaralLab/yucheng/HMM_histone
#     bash run_all.sh
#
# It will:
#   1. create the 'hmmchip' conda env if it doesn't exist yet
#   2. download the ENCODE data (runs right here; a few minutes)
#   3. show you the chosen ENCODE accessions and ask you to confirm
#   4. submit 01 -> 02 -> 03 as a SLURM dependency chain
#      (each step only starts if the previous one succeeded)
#
# After that you're done until the email says the jobs finished.
# ============================================================

set -euo pipefail

base_dir="/projects/b1042/AmaralLab/yucheng/HMM_histone"
cd "$base_dir"

module load python-miniconda3

# --- 1. conda env (one-time, auto-skipped if present) ---
if conda env list | awk '{print $1}' | grep -qx hmmchip; then
    echo "[env] 'hmmchip' env already exists, skipping setup."
else
    echo "[env] Creating 'hmmchip' conda env (one-time, ~5 min)..."
    conda create -n hmmchip python=3.9 numpy pandas matplotlib seaborn -y
    conda run -n hmmchip pip install pyhhmm pyBigWig
    echo "[env] Done."
fi

# --- 2. download data (runs here, not via sbatch) ---
manifest="data/ENCSR000DWD.manifest.tsv"
if [ -s "$manifest" ] && ls data/K562_H3K4me3/bams/*.bam >/dev/null 2>&1; then
    echo "[download] Data already present, skipping download."
else
    echo "[download] Downloading ENCODE data..."
    bash scripts/00_download_data.sh
fi

# --- 3. manifest checkpoint ---
echo
echo "================ ENCODE files selected ================"
column -t -s $'\t' "$manifest"
echo "========================================================"
echo "Sanity check: 2 BAMs (output_type 'alignments' = the FILTERED"
echo "alignments, one per isogenic replicate), 1 fold-change bigWig,"
echo "1 narrowPeak -- all GRCh38."
echo
read -r -p "Look right? Submit the SLURM jobs? [y/N] " ans
case "$ans" in
    [yY]*) ;;
    *) echo "Stopped. Nothing was submitted."; exit 0 ;;
esac

# --- 4. submit the dependency chain ---
jid1=$(sbatch --parsable scripts/01_raw_coverage_bw.sh)
echo "[submit] 01_raw_coverage_bw.sh      -> job $jid1"

jid2=$(sbatch --parsable --dependency=afterok:"$jid1" scripts/02_pyhhmm.sh)
echo "[submit] 02_pyhhmm.sh (after $jid1) -> job $jid2"

jid3=$(sbatch --parsable --dependency=afterok:"$jid2" scripts/03_aic_bic.sh)
echo "[submit] 03_aic_bic.sh (after $jid2) -> job $jid3"

echo
echo "All queued. Check progress anytime with:"
echo "    squeue -u \$USER"
echo "You'll get an email at each job's END or FAIL."
echo
echo "When job $jid3 finishes, results are in:"
echo "    data/pyhhmm/K562_H3K4me3_raw/"
echo "    - aic_bic_table.chr14.csv + aic_bic_barplots.chr14.png  (pick best k)"
echo "    - *<k>states_regions.bed                                (load in IGV)"
