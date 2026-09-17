#!/bin/bash
#SBATCH -A b1042
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --time=48:00:00
#SBATCH --mem=10G
#SBATCH --partition  genomics
#SBATCH --job-name=raw_coverage_1bp_and_100bp_mean_H3K4me3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yuchengjin2027@u.northwestern.edu

# Adapted from Liebe's raw_coverage_bw.sh .
# Produces the raw total coverage bigWigs that feed the Gaussian HMM:
#   1) bamCoverage at 1 bp, no normalization
#   2) bigwigCompare mean-binned to 100 bp
# This mirrors Liebe's method 2 for 100 bp bins ("binned outputs of 1 bp
# separately to get average pileup per bin", handover PDF p.3).
#
# ChIP-seq-specific notes vs. the ATAC original (parameters kept identical
# on purpose):
#   --extendReads 150      : ENCSR000DWD is single-end (Illumina GA, short
#                            reads); SE reads must be extended toward the
#                            fragment size. H3K4me3 sonication fragments are
#                            ~150-250 bp, so 150 remains reasonable.
#   --ignoreDuplicates and --minMappingQuality 30 : the ENCODE "alignments"
#                            BAMs are already filtered (dedup + MAPQ), so
#                            these are mostly no-ops here -- kept to match
#                            Liebe's script exactly.

set -euo pipefail

# NOTE: no conda env here on purpose. This script only uses module-provided
# tools (bamCoverage, bigwigCompare, samtools). Activating hmmchip would put
# its Python first on PATH and could shadow the deeptools module's
# interpreter -- same interpreter-shadowing failure mode as on the figure
# reproduction track. The hmmchip env is needed only in scripts 02 and 03.
module load deeptools
module load samtools


#---paths---
base_dir="/projects/b1042/AmaralLab/yucheng/HMM_histone"
in_dir="${base_dir}/data/K562_H3K4me3/bams"
out_dir="${base_dir}/data/K562_H3K4me3/raw_coverage_bw"

mkdir -p "$out_dir"


#---run bamCoverage---
binsize=1

for bam in "$in_dir"/*.bam; do
    [ -e "$bam" ] || continue   # skip if no .bam files
    base=$(basename "$bam" .bam)
    out_bw="$out_dir/${base}.rawCoverage.${binsize}bp.bw"

    echo "Processing $bam"
    echo "Output: $out_bw"

    if [ ! -f "${bam}.bai" ]; then
    	echo "Indexing $bam"
    	samtools index "$bam"
    fi

    bamCoverage \
        -b "$bam" \
        -o "$out_bw" \
        --binSize "$binsize" \
        --normalizeUsing None \
        --extendReads 150 \
        --ignoreDuplicates \
        --minMappingQuality 30 \
        -p 10

    echo "Done: $base"
    echo
done


#---get averaged bin by 100bp---
bin100=100

for bw in "$out_dir"/*.rawCoverage.${binsize}bp.bw; do
    [ -e "$bw" ] || continue

    base=$(basename "$bw" .rawCoverage.1bp.bw)
    out_bw_100bp="$out_dir/${base}.rawCoverage.${bin100}bp.mean.bw"

    echo "Binning $bw -> $out_bw_100bp"

    bigwigCompare \
        -b1 "$bw" \
        -b2 "$bw" \
        --operation mean \
        --binSize "$bin100" \
        --skipNonCoveredRegions \
        -o "$out_bw_100bp"

    echo "Done binning: $base"
    echo
done


echo "All done."
