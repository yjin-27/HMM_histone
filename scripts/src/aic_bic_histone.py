"""
AIC/BIC model selection across the 2-7 state HMMs.

Reproduces the "AIC BIC num states" section of Liebe's pyhhmm.ipynb as a
batch scripts. Parses the final logL
out of each *parameters.txt written by gaussian_hmm_histone.py, computes
AIC and BIC per (replicate, k), and writes:
    aic_bic_table.csv
    aic_bic_barplots.png   (per-replicate + averaged, matching the notebook)

Formulas match the notebook exactly:
    n_params = (k - 1)            initial probs
             + k * (k - 1)        transition matrix
             + k * n_emissions    means
             + k * n_emissions    diagonal variances
    AIC = 2 * n_params - 2 * logL
    BIC = n_params * ln(n) - 2 * logL
where n = number of 100 bp bins on the analyzed chromosome.
"""

import glob
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import seaborn as sns

base_dir = '/projects/b1042/AmaralLab/yucheng/HMM_histone'
hmm_dir = os.path.join(base_dir, 'data/pyhhmm/K562_H3K4me3_raw')
bw_dir = os.path.join(base_dir, 'data/K562_H3K4me3/raw_coverage_bw')

resolution = 100
chrom = os.environ.get('CHROM', 'chr14')
n_emissions = 1


def calculate_n_params(k, n_emissions=1, covariance_type='diagonal'):
    assert covariance_type == 'diagonal'
    return (k - 1) + k * (k - 1) + k * n_emissions + k * n_emissions


def get_sample_size(bw_path, chrom, resolution):
    bw = pyBigWig.open(bw_path)
    if chrom not in bw.chroms():
        bw.close()
        return None
    n = int(np.ceil(bw.chroms()[chrom] / resolution))
    bw.close()
    return n


def final_logL(params_file):
    """Last 'Iteration i: <logL>' line in the parameters file."""
    ll = None
    with open(params_file) as f:
        for line in f:
            m = re.match(r'\s*Iteration \d+: (-?\d+\.?\d*)', line)
            if m:
                ll = float(m.group(1))
    return ll


# --- sample size n from any replicate's bigWig ---
ref_bws = sorted(glob.glob(os.path.join(bw_dir, '*.rawCoverage.100bp.mean.bw')))
assert ref_bws, f"no raw coverage bigWigs in {bw_dir}"
n_samples = get_sample_size(ref_bws[0], chrom, resolution)
assert n_samples, f"{chrom} not found in {ref_bws[0]}"
print(f"n (bins on {chrom} at {resolution}bp): {n_samples}")

# --- collect logL per (replicate, k) ---
rows = []
pattern = os.path.join(hmm_dir, f'*.{chrom}.*states_parameters.txt')
for pf in sorted(glob.glob(pattern)):
    base = os.path.basename(pf)
    m = re.match(r'.+\.(ENCFF\w+)\.(\d+)states_parameters\.txt', base)
    if not m:
        print(f"skipping unparsable filename: {base}")
        continue
    rep, k = m.group(1), int(m.group(2))
    ll = final_logL(pf)
    if ll is None:
        print(f"no logL found in {base}, skipping")
        continue
    rows.append({'replicate': rep, 'k': k, 'logL': ll})

df = pd.DataFrame(rows).sort_values(['replicate', 'k'])
assert not df.empty, f"no parameter files matched {pattern}"

# --- AIC / BIC ---
df['n_params'] = df['k'].apply(lambda k: calculate_n_params(k, n_emissions, 'diagonal'))
df['AIC'] = 2 * df['n_params'] - 2 * df['logL']
df['BIC'] = df['n_params'] * np.log(n_samples) - 2 * df['logL']

out_csv = os.path.join(hmm_dir, f'aic_bic_table.{chrom}.csv')
df.to_csv(out_csv, index=False)
print(df.to_string(index=False))
print(f"\nWrote {out_csv}")

# --- summary across replicates ---
summary = df.groupby('k')[['AIC', 'BIC']].agg(['mean', 'std']).reset_index()
summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]

# --- plots: one panel per replicate + averaged panel ---
reps = df['replicate'].unique()
n_panels = len(reps) + 1
fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6), squeeze=False)
axes = axes[0]

for ax, rep in zip(axes, reps):
    rep_df = df[df['replicate'] == rep].sort_values('k')
    melted = pd.melt(rep_df, id_vars=['k'], value_vars=['AIC', 'BIC'],
                     var_name='Metric', value_name='Value')
    sns.barplot(x='k', y='Value', hue='Metric', data=melted, ax=ax)
    ax.set_title(f'AIC & BIC vs. k -- {rep} ({chrom})')
    ax.set_xlabel('Number of States (k)')
    ax.set_ylabel('Information Criterion Value')

ax = axes[-1]
melted = pd.melt(summary, id_vars=['k'], value_vars=['AIC_mean', 'BIC_mean'],
                 var_name='Metric', value_name='Value')
sns.barplot(x='k', y='Value', hue='Metric', data=melted, ax=ax)
ax.set_title(f'Average AIC & BIC across replicates ({chrom})')
ax.set_xlabel('Number of States (k)')
ax.set_ylabel('Mean Information Criterion Value')

plt.tight_layout()
out_png = os.path.join(hmm_dir, f'aic_bic_barplots.{chrom}.png')
plt.savefig(out_png, dpi=150)
print(f"Wrote {out_png}")

best_aic = df.loc[df.groupby('replicate')['AIC'].idxmin(), ['replicate', 'k', 'AIC']]
best_bic = df.loc[df.groupby('replicate')['BIC'].idxmin(), ['replicate', 'k', 'BIC']]
print("\nBest k by AIC per replicate:\n", best_aic.to_string(index=False))
print("\nBest k by BIC per replicate:\n", best_bic.to_string(index=False))
