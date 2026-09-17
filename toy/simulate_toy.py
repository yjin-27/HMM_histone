#!/usr/bin/env python
"""
Stage 1 toy dataset for the H3K4me3 Gaussian HMM.

Simulates a 3-state Gaussian HMM with sticky transitions, writes the result as a
bigWig in the same layout the real pipeline consumes, and saves ground truth
(state path, generating parameters, exact log-likelihood) alongside it.

The pipeline applies log1p to whatever it reads out of the bigWig, so emissions
are generated in log1p space and inverted with expm1 before writing. This means
the HMM sees exactly the means/sigmas specified here.

Variants (see --variant):
    A  means (2, 4, 8), sigma 1.00   -> sigma > 1/sqrt(2*pi), expect POSITIVE AIC/BIC
    B  means (2, 4, 8), sigma 0.05   -> sigma < 1/sqrt(2*pi), expect NEGATIVE AIC/BIC
    C  means (0, 0.5, 3), sigma 0.05 -> realistic closed/mid/open

A and B differ only in sigma, so they isolate the sign question.

Usage:
    python simulate_toy.py --variant A --outdir data/toy
    python simulate_toy.py --variant B --outdir data/toy
    python simulate_toy.py --variant C --outdir data/toy --n-bins 1000
    python simulate_toy.py --variant C --outdir data/toy --zero-frac 0.6
"""

import argparse
import json
import os

import numpy as np

# sigma below this makes the Gaussian density exceed 1 at the mean,
# which makes log-likelihood contributions positive and AIC/BIC negative.
SIGMA_DENSITY_THRESHOLD = 1.0 / np.sqrt(2.0 * np.pi)  # ~0.3989

VARIANTS = {
    'A': {'means': [2.0, 4.0, 8.0], 'sigmas': [1.00, 1.00, 1.00]},
    'B': {'means': [2.0, 4.0, 8.0], 'sigmas': [0.05, 0.05, 0.05]},
    'C': {'means': [0.0, 0.5, 3.0], 'sigmas': [0.05, 0.05, 0.05]},
}


def build_transmat(k, stay_prob):
    """Sticky transition matrix: stay_prob on the diagonal, rest spread evenly."""
    off = (1.0 - stay_prob) / (k - 1)
    T = np.full((k, k), off)
    np.fill_diagonal(T, stay_prob)
    return T


def stationary_distribution(transmat):
    """Left eigenvector of the transition matrix with eigenvalue 1."""
    vals, vecs = np.linalg.eig(transmat.T)
    idx = np.argmin(np.abs(vals - 1.0))
    pi = np.real(vecs[:, idx])
    return pi / pi.sum()


def simulate(means, sigmas, transmat, n_bins, seed, zero_frac=0.0):
    """
    Draw a state path from the Markov chain, then emit a Gaussian per bin.

    Returns (emissions_log1p, true_states). Emissions are clipped at 0 because
    coverage cannot be negative -- the same floor the real data has.
    """
    rng = np.random.default_rng(seed)
    k = len(means)
    pi = stationary_distribution(transmat)

    states = np.empty(n_bins, dtype=int)
    states[0] = rng.choice(k, p=pi)
    for t in range(1, n_bins):
        states[t] = rng.choice(k, p=transmat[states[t - 1]])

    emissions = rng.normal(loc=np.asarray(means)[states],
                           scale=np.asarray(sigmas)[states])
    emissions = np.clip(emissions, 0.0, None)

    if zero_frac > 0:
        # Optional Stage-4 probe: force a fraction of bins to exact zero,
        # mimicking unmappable / uncovered regions in real coverage tracks.
        mask = rng.random(n_bins) < zero_frac
        emissions[mask] = 0.0

    return emissions, states


def forward_loglik(x, means, sigmas, transmat, pi):
    """
    Exact marginal log-likelihood under the GIVEN parameters, via the scaled
    forward algorithm. This is the reference value that a correctly fitted
    model should slightly exceed. If Baum-Welch lands below this, it found a
    bad local optimum.
    """
    n, k = len(x), len(means)
    means = np.asarray(means)
    sigmas = np.asarray(sigmas)

    # log emission densities, shape (n, k)
    z = (x[:, None] - means[None, :]) / sigmas[None, :]
    log_b = -0.5 * z ** 2 - np.log(sigmas)[None, :] - 0.5 * np.log(2 * np.pi)
    b = np.exp(log_b)

    loglik = 0.0
    alpha = pi * b[0]
    c = alpha.sum()
    loglik += np.log(c)
    alpha /= c
    for t in range(1, n):
        alpha = (alpha @ transmat) * b[t]
        c = alpha.sum()
        loglik += np.log(c)
        alpha /= c
    return loglik


def write_bigwig(path, chrom, values, resolution):
    """Write one value per bin as fixed-width intervals, in RAW (expm1) space."""
    import pyBigWig  # imported here so the rest of the script runs without it

    raw = np.expm1(values).astype(np.float64)
    raw = np.clip(raw, 0.0, None)
    chrom_len = len(values) * resolution

    bw = pyBigWig.open(path, 'w')
    bw.addHeader([(chrom, chrom_len)])
    starts = np.arange(len(values), dtype=np.int64) * resolution
    bw.addEntries([chrom] * len(values),
                  starts.tolist(),
                  ends=(starts + resolution).tolist(),
                  values=raw.tolist())
    bw.close()
    return chrom_len


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--variant', choices=sorted(VARIANTS), default='C')
    p.add_argument('--n-bins', type=int, default=1000)
    p.add_argument('--resolution', type=int, default=100)
    p.add_argument('--chrom', default='chrTOY')
    p.add_argument('--stay-prob', type=float, default=0.95,
                   help='diagonal of the transition matrix; 0.95 -> mean run ~20 bins')
    p.add_argument('--zero-frac', type=float, default=0.0,
                   help='fraction of bins forced to exact zero (Stage 4 probe)')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--outdir', default='.')
    p.add_argument('--tag', default=None, help='override the output basename')
    p.add_argument('--no-bigwig', action='store_true',
                   help='skip bigWig output (useful if pyBigWig is unavailable)')
    args = p.parse_args()

    spec = VARIANTS[args.variant]
    means, sigmas = spec['means'], spec['sigmas']
    k = len(means)

    transmat = build_transmat(k, args.stay_prob)
    pi = stationary_distribution(transmat)

    emissions, states = simulate(means, sigmas, transmat, args.n_bins,
                                 args.seed, zero_frac=args.zero_frac)

    os.makedirs(args.outdir, exist_ok=True)
    tag = args.tag or (f'toy{args.variant}_T{args.n_bins}_seed{args.seed}'
                       + (f'_z{args.zero_frac:g}' if args.zero_frac else ''))

    # ground truth -- everything Stage 2 needs to score the fit
    npz_path = os.path.join(args.outdir, f'{tag}.truth.npz')
    np.savez(npz_path,
             emissions_log1p=emissions,
             true_states=states,
             means=np.asarray(means),
             sigmas=np.asarray(sigmas),
             transmat=transmat,
             startprob=pi)

    true_ll = forward_loglik(emissions, means, sigmas, transmat, pi)
    per_bin = true_ll / args.n_bins

    meta = {
        'variant': args.variant, 'chrom': args.chrom,
        'n_bins': args.n_bins, 'resolution': args.resolution,
        'k_true': k, 'means': means, 'sigmas': sigmas,
        'stay_prob': args.stay_prob, 'zero_frac': args.zero_frac,
        'seed': args.seed,
        'true_loglik': true_ll, 'true_loglik_per_bin': per_bin,
        'sigma_density_threshold': SIGMA_DENSITY_THRESHOLD,
        'expect_negative_aic': bool(per_bin > 0),
    }

    bw_path = None
    if not args.no_bigwig:
        bw_path = os.path.join(args.outdir, f'{tag}.rawCoverage.{args.resolution}bp.mean.bw')
        meta['chrom_length'] = write_bigwig(bw_path, args.chrom, emissions, args.resolution)

    with open(os.path.join(args.outdir, f'{tag}.meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    # ---- report ----
    occ = np.bincount(states, minlength=k) / args.n_bins
    print(f'variant {args.variant}  T={args.n_bins}  seed={args.seed}')
    print(f'  means           {means}')
    print(f'  sigmas          {sigmas}   (density>1 threshold: {SIGMA_DENSITY_THRESHOLD:.4f})')
    print(f'  occupancy       {np.round(occ, 3).tolist()}')
    print(f'  emission range  [{emissions.min():.3f}, {emissions.max():.3f}] (log1p space)')
    print(f'  true logL       {true_ll:,.1f}   ({per_bin:+.4f} per bin)')
    print(f'  PREDICTION      AIC/BIC will be '
          f'{"NEGATIVE" if per_bin > 0 else "POSITIVE"} at the true parameters')
    print(f'  truth           {npz_path}')
    if bw_path:
        print(f'  bigWig          {bw_path}')


if __name__ == '__main__':
    main()
