# HMM segmentation of histone ChIP-seq (H3K4me3) — Quest workflow

Testing Liebe's Gaussian-HMM raw-coverage segmentation (from the ATAC-seq
project) on a histone ChIP-seq dataset, per discussion with Maalavika.

**Dataset:** [ENCSR000DWD](https://www.encodeproject.org/experiments/ENCSR000DWD/) —
H3K4me3 Histone ChIP-seq, K562, two **isogenic replicates** (Stamatoyannopoulos
lab, ENCODE2). Control experiment (used by ENCODE's fold-change track, *not* by
our HMM): ENCSR000DWA.

**Question:** does the HMM segmentation of *raw coverage* (no control, no
enrichment model) recover the same regions as ENCODE's control-normalized
fold-change bigWig + narrowPeak calls? And do the two isogenic replicates give
concordant segmentations?

---

## Directory layout (on Quest)

```
/projects/b1042/AmaralLab/yucheng/HMM_histone/
├── README.md
├── data/
│   ├── hg38.chrom.sizes
│   ├── ENCSR000DWD.manifest.tsv          # which ENCODE files were chosen
│   ├── K562_H3K4me3/
│   │   ├── bams/                         # 2 filtered alignment BAMs
│   │   ├── raw_coverage_bw/              # 1bp + 100bp-mean bigWigs
│   │   └── encode_reference/             # ENCODE FC bigWig + narrowPeak
│   └── pyhhmm/K562_H3K4me3_raw/          # HMM outputs
└── scripts/
    ├── 00_download_data.sh
    ├── 01_raw_coverage_bw.sh
    ├── 02_pyhhmm.sh
    ├── 03_aic_bic.sh
    └── src/
        ├── gaussian_hmm_histone.py
        └── aic_bic_histone.py
```

## 0. One-time environment setup

Liebe's `jupyter-kernel-py37` env is hers; make your own:

```bash
module load python-miniconda3
conda create -n hmmchip python=3.9 numpy pandas matplotlib seaborn -y
conda activate hmmchip
pip install pyhhmm pyBigWig
```

(`samtools` comes with the deeptools module load in script 01, matching
Liebe's setup. Edit the `mail-user` placeholder in each script.)

## 1. Download data — `00_download_data.sh`

Queries the ENCODE REST API and picks, for GRCh38:

- the **two filtered alignment BAMs** (`output_type = "alignments"` — this is
  ENCODE's name for the *filtered* BAMs; the unfiltered ones are
  `"unfiltered alignments"`), one per isogenic replicate;
- the **fold change over control** bigWig (pooled-replicate track preferred);
- a **narrowPeak** file (replicated peaks preferred — the histone pipeline's
  equivalent of the IDR peaks mentioned in the handover for TF/ATAC data).

Every chosen accession is echoed and written to
`data/ENCSR000DWD.manifest.tsv` — **check it against the portal file table
before proceeding** (File format = bam, Output type = alignments, as in the
handover PDF's data-selection instructions).

Run inside tmux on a login node, or `sbatch` it.

## 2. Raw coverage bigWigs — `01_raw_coverage_bw.sh`

Direct adaptation of Liebe's `raw_coverage_bw.sh`: `bamCoverage` at 1 bp with
`--normalizeUsing None --extendReads 150 --ignoreDuplicates
--minMappingQuality 30`, then `bigwigCompare --operation mean --binSize 100`
(her "method 2" for 100 bp bins). Parameters intentionally unchanged so the
histone results are directly comparable to her ATAC runs. Notes:

- ENCSR000DWD is **single-end**, so `--extendReads 150` extends reads toward
  fragment size (appropriate for ~150–250 bp sonication fragments).
- The ENCODE filtered BAMs are already dedup/MAPQ-filtered, so the
  `--ignoreDuplicates` / `-q 30` flags are mostly redundant here — kept for
  consistency.

`sbatch scripts/01_raw_coverage_bw.sh` (genome-wide 1 bp coverage: give it the
full job, not a login node).

## 3. Train HMMs — `02_pyhhmm.sh`

Runs `src/gaussian_hmm_histone.py`, adapted from `gaussian_hmm_raw.py`:

- both replicates processed in parallel (multiprocessing, one process each);
- **chr14 by default** (Liebe's main chromosome; `export CHROM=chr20` etc. to
  override — H3K4me3 is promoter-centric, so any gene-dense chromosome works);
- observations = 100 bp binned mean coverage, `log1p`-transformed;
- sweeps **k = 2…7 states** (restored from the handover; Liebe's checked-in
  copy had been narrowed to k = 7 only);
- per (replicate, k) writes `*parameters.txt`, `*regions.bed` (IGV,
  color-coded low→high signal: purple→crimson), `*transitions.bed`.

## 4. Model selection — `03_aic_bic.sh`

Batch-script port of the notebook's "AIC BIC num states" section: parses final
logL from each parameters file, computes AIC/BIC with the same n_params
formula, writes `aic_bic_table.<chrom>.csv` + `aic_bic_barplots.<chrom>.png`,
and prints the best k per replicate. The paper-goal check: is the optimal k
still > 2 for a histone mark?

## 5. Visual comparison in IGV

Load in IGV (hg38 genome), for the chromosome you segmented:

1. `raw_coverage_bw/*.rawCoverage.100bp.mean.bw` (your HMM input, per rep)
2. `encode_reference/*.fc.signal.bigwig` (ENCODE fold change over control)
3. `encode_reference/*.narrowPeak.bed.gz` (ENCODE peak calls)
4. `pyhhmm/K562_H3K4me3_raw/*<best_k>states_regions.bed` for **both**
   replicates

What to look for:

- **HMM vs. ENCODE:** do the high-emission states (red/crimson) coincide with
  narrowPeaks and FC-track summits? H3K4me3 is sharp and promoter-focused, so
  expect tight high-state blocks at TSSs. Where they disagree, is it the
  control normalization (present in FC, absent in raw coverage) doing the work?
- **Rep1 vs. Rep2:** isogenic replicates → segmentations should be highly
  concordant. Discordance is a reproducibility signal about the method, not
  the biology.
- **Intermediate states:** does a histone mark also support >2 meaningful
  states (e.g., broad weak H3K4me3 shoulders vs. sharp summits), mirroring the
  ATAC conclusion?

Optional quantification (same idea as eyeballing, but numeric):

```bash
module load bedtools/2.29.2
# concordance of top-state calls between reps, and vs ENCODE peaks
grep "state_<top>" rep1_regions.bed | bedtools merge -i - > rep1_top.bed
grep "state_<top>" rep2_regions.bed | bedtools merge -i - > rep2_top.bed
bedtools jaccard -a rep1_top.bed -b rep2_top.bed
zcat encode_reference/*.narrowPeak.bed.gz | awk -v c=$CHROM '$1==c' | \
    sort -k1,1 -k2,2n | bedtools jaccard -a rep1_top.bed -b -
```

## Quickstart (the only commands you need)

```bash
cd /projects/b1042/AmaralLab/yucheng
# put the HMM_histone folder here, then:
cd HMM_histone
tmux new -s hmm          # optional but recommended
bash run_all.sh
```

`run_all.sh` creates the conda env if needed, downloads the data, shows you
the chosen ENCODE accessions for a quick y/N confirmation, then submits
01 → 02 → 03 as a SLURM dependency chain. You'll get an email when each job
ends. Everything below documents what the individual steps do; you don't need
to run them by hand.

When the last job finishes: pick the best k from
`data/pyhhmm/K562_H3K4me3_raw/aic_bic_barplots.chr14.png`, rsync the bigWigs +
BEDs to your laptop, and do the IGV comparison in section 5.
