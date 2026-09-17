# HMM segmentation of histone ChIP-seq coverage

Gaussian HMM chromatin-state segmentation applied to H3K4me3 ChIP-seq in K562
(ENCODE [ENCSR000DWD](https://www.encodeproject.org/experiments/ENCSR000DWD/),
two isogenic replicates). Adapted from Liebe's ATAC-seq raw-coverage pipeline.

Peak callers such as MACS2 were designed for transcription factors and treat
enrichment as binary. Histone marks are domain-structured with continuous
signal intensity, so an HMM — which segments the whole genome and models
spatial autocorrelation through its transition matrix — is a closer structural
match. This repository tests whether the ATAC pipeline transfers, and includes
a synthetic ground-truth framework for validating the fitting procedure.

## Files

### `scripts/`

| file | what it does |
|---|---|
| `00_download_data.sh` | Queries the ENCODE REST API for the two filtered alignment BAMs, the fold-change bigWig, and a narrowPeak file. Writes chosen accessions to a manifest for checking against the portal. |
| `01_raw_coverage_bw.sh` | `bamCoverage` at 1 bp (`--normalizeUsing None --extendReads 150 --ignoreDuplicates --minMappingQuality 30`), then `bigwigCompare --operation mean --binSize 100`. Parameters unchanged from Liebe's so results stay comparable to the ATAC runs. |
| `02_pyhhmm.sh` | SLURM wrapper for the HMM training script. Resolves the conda env's interpreter explicitly, since `module load python-miniconda3` can shadow it. |
| `03_aic_bic.sh` | Parses final log-likelihoods from the parameter files, computes AIC/BIC, writes a table and bar plots. |
| `src/gaussian_hmm_histone.py` | The pipeline: reads binned coverage from bigWig, `log1p`-transforms, sweeps k = 2-7 states, writes parameters + IGV-colored state BEDs + transition BEDs. |
| `src/aic_bic_histone.py` | Model-selection computation called by `03`. |

### `toy/`

Synthetic data with known ground truth, used to separate genuine model
behaviour from artifacts of the optimizer.

| file | what it does |
|---|---|
| `simulate_toy.py` | Generates a 3-state Gaussian HMM with sticky transitions and writes it as a bigWig, so toy data enters the pipeline at the same point real data does. Saves the true state path, generating parameters, and the exact forward-algorithm log-likelihood at those parameters. |
| `score_toy.py` | Reads the fitted parameter files, compares against ground truth, and reports AIC/BIC, parameter recovery, and any place where log-likelihood decreased with k. |

Three variants: **A** and **B** share means (2, 4, 8) and differ only in
emission sigma (1.0 vs 0.05), isolating the AIC/BIC sign question. **C** uses
(0, 0.5, 3.0), closer to closed / intermediate / open chromatin.

## Configuration

`gaussian_hmm_histone.py` reads environment variables (defaults in parentheses):

| var | meaning | default |
|---|---|---|
| `DATA_DIR` | input bigWig directory | `data/K562_H3K4me3/raw_coverage_bw` |
| `OUTPUT_DIR` | where results are written | `data/pyhhmm/K562_H3K4me3_raw` |
| `CHROM` | chromosome to segment | `chr14` |
| `CELL_LINE`, `MARK` | output filename prefixes | `K562`, `H3K4me3` |
| `N_INIT` | random restarts per model | `1` |
| `N_ITER` | max EM iterations per restart | `20` |
| `NPROC` | worker processes | all cores |

## Running the toy validation

    python toy/simulate_toy.py --variant C --seed 0 --chrom chrTOY --outdir toy/data

    DATA_DIR=$PWD/toy/data OUTPUT_DIR=$PWD/toy/pyhhmm CHROM=chrTOY \
      CELL_LINE=TOY MARK=sim N_INIT=10 N_ITER=200 \
      python scripts/src/gaussian_hmm_histone.py

    python toy/score_toy.py --params-dir toy/pyhhmm --truth-dir toy/data


## Environment

    module load python-miniconda3
    conda create -n hmmchip python=3.9 numpy pandas matplotlib seaborn -y
    conda activate hmmchip
    pip install pyhhmm pyBigWig

Data directories are gitignored. Toy data regenerates exactly from `--seed`.
