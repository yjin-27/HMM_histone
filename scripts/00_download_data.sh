#!/bin/bash
# Hard-coded downloader for the H3K4me3 HMM test (ENCSR000DWD).
# ----------------------------------------------------------------
# EDIT THE FOUR ACCESSIONS BELOW, then run:   bash scripts/00_download_data.sh
#
# Read them off https://www.encodeproject.org/experiments/ENCSR000DWD/
# file table (default view = default analysis), GRCh38:
#   - the two "alignments" BAMs (Rep 1 and Rep 2)
#   - the pooled "signal p-value" bigWig (replicates "1, 2")
#   - the "replicated peaks" bed narrowPeak
# ----------------------------------------------------------------

REP1_BAM="ENCFF185YRK"      # Rep 1 filtered alignments (bam, GRCh38)
REP2_BAM="ENCFF955AMI"      # Rep 2 filtered alignments (bam, GRCh38)
FC_BIGWIG="ENCFF012YLR"     # signal p-value bigWig (pooled, GRCh38, portal default)
NARROWPEAK="ENCFF148POZ"    # replicated peaks (bed narrowPeak, pooled, GRCh38)

set -euo pipefail

base_dir="/projects/b1042/AmaralLab/yucheng/HMM_histone"
data_dir="${base_dir}/data"
bam_dir="${data_dir}/K562_H3K4me3/bams"
encode_dir="${data_dir}/K562_H3K4me3/encode_reference"
mkdir -p "$bam_dir" "$encode_dir"

for acc in "$REP1_BAM" "$REP2_BAM" "$FC_BIGWIG" "$NARROWPEAK"; do
    if [[ "$acc" == ENCFFXXXXXX ]]; then
        echo "ERROR: edit the accessions at the top of this script first."
        exit 1
    fi
done

dl () {  # dl <accession> <extension> <output_path>
    local acc=$1 ext=$2 out=$3
    if [ -s "$out" ]; then
        echo "Already present: $out"
    else
        echo "Downloading $acc -> $out"
        curl -sSL "https://www.encodeproject.org/files/${acc}/@@download/${acc}.${ext}" -o "$out"
    fi
}

dl "$REP1_BAM"   bam            "$bam_dir/${REP1_BAM}.bam"
dl "$REP2_BAM"   bam            "$bam_dir/${REP2_BAM}.bam"
dl "$FC_BIGWIG"  bigWig         "$encode_dir/${FC_BIGWIG}.pval.signal.bigwig"
dl "$NARROWPEAK" bed.gz         "$encode_dir/${NARROWPEAK}.narrowPeak.bed.gz"

# chrom sizes
chrsz="${data_dir}/hg38.chrom.sizes"
[ -s "$chrsz" ] || curl -sSL \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" -o "$chrsz"

# manifest (kept so run_all.sh's checkpoint still works)
manifest="${data_dir}/ENCSR000DWD.manifest.tsv"
{
    echo -e "role\taccession\toutput_type\tassembly\treplicates\thref"
    echo -e "bam\t${REP1_BAM}\talignments\tGRCh38\t1\thttps://www.encodeproject.org/files/${REP1_BAM}/"
    echo -e "bam\t${REP2_BAM}\talignments\tGRCh38\t2\thttps://www.encodeproject.org/files/${REP2_BAM}/"
    echo -e "fc_bigwig\t${FC_BIGWIG}\tsignal p-value\tGRCh38\t1,2\thttps://www.encodeproject.org/files/${FC_BIGWIG}/"
    echo -e "narrowpeak\t${NARROWPEAK}\treplicated peaks\tGRCh38\t1,2\thttps://www.encodeproject.org/files/${NARROWPEAK}/"
} > "$manifest"

echo
echo "Done. Files:"
ls -lh "$bam_dir" "$encode_dir"
echo "Manifest: $manifest"
