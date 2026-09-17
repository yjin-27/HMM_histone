#!/usr/bin/env python
"""
Score the toy HMM fits against ground truth.

Reads the *states_parameters.txt files written by gaussian_hmm_histone.py,
pairs each with its .meta.json / .truth.npz, and reports:
  - fitted logL vs. the exact logL at the TRUE parameters
  - AIC / BIC and which k they select
  - whether logL is monotone non-decreasing in k (it must be, in theory)
  - recovered means/variances at k=3 vs. the generating values

Usage:
    python score_toy.py --params-dir toy/pyhhmm --truth-dir toy/data
"""

import argparse
import glob
import json
import os
import re

import numpy as np


def parse_params(path):
    """Pull k, final logL, means and variances out of one parameter file."""
    fname = os.path.basename(path)

    m = re.search(r'\.(\d+)states_parameters\.txt$', fname)
    if not m:
        return None
    k = int(m.group(1))

    # sample name is the field before '.<k>states_parameters.txt'
    sample = fname[:m.start()].split('.')[-1]

    text = open(path).read()

    it_lines = [l for l in text.splitlines() if l.strip().startswith('Iteration')]
    if not it_lines:
        return None
    m_ll = re.search(r':\s*(-?[\d.]+)', it_lines[-1])
    if not m_ll:
        return None

    return {
        'sample': sample,
        'k': k,
        'logL': float(m_ll.group(1)),
        'n_iter_run': len(it_lines),
        'means': [float(v) for v in re.findall(r'Mean signal = (-?[\d.]+)', text)],
        'variances': [float(v) for v in re.findall(r'Variance = (-?[\d.]+)', text)],
    }


def n_params(k, n_emissions=1):
    """(k-1) initial + k(k-1) transitions + k means + k diagonal variances."""
    return (k - 1) + k * (k - 1) + k * n_emissions + k * n_emissions


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--params-dir', default='toy/pyhhmm')
    p.add_argument('--truth-dir', default='toy/data')
    args = p.parse_args()

    fits = [r for r in (parse_params(f) for f in
                        sorted(glob.glob(os.path.join(args.params_dir,
                                                      '*states_parameters.txt'))))
            if r]
    if not fits:
        raise SystemExit(f'no parameter files found in {args.params_dir}')

    by_sample = {}
    for r in fits:
        by_sample.setdefault(r['sample'], []).append(r)

    for sample in sorted(by_sample):
        rows = sorted(by_sample[sample], key=lambda r: r['k'])

        meta_path = os.path.join(args.truth_dir, f'{sample}.meta.json')
        meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
        true_ll = meta.get('true_loglik')
        n_bins = meta.get('n_bins', 1000)
        k_true = meta.get('k_true')

        print(f'\n{"=" * 78}')
        print(f'{sample}   T={n_bins}  k_true={k_true}  '
              f'means={meta.get("means")}  sigmas={meta.get("sigmas")}')
        if true_ll is not None:
            print(f'  exact logL at TRUE parameters: {true_ll:>12,.2f}   '
                  f'({true_ll / n_bins:+.4f} / bin)')
        print(f'{"=" * 78}')

        hdr = f'{"k":>3} {"logL":>12} {"/bin":>9} {"vs true":>11} {"p":>4} {"AIC":>12} {"BIC":>12} {"iters":>6}'
        print(hdr)
        print('-' * len(hdr))

        best_aic = min(rows, key=lambda r: 2 * n_params(r['k']) - 2 * r['logL'])
        best_bic = min(rows, key=lambda r: n_params(r['k']) * np.log(n_bins) - 2 * r['logL'])

        for r in rows:
            pk = n_params(r['k'])
            aic = 2 * pk - 2 * r['logL']
            bic = pk * np.log(n_bins) - 2 * r['logL']
            delta = f'{r["logL"] - true_ll:+,.1f}' if true_ll is not None else '--'
            flag = ''
            if r is best_aic:
                flag += ' <AIC'
            if r is best_bic:
                flag += ' <BIC'
            print(f'{r["k"]:>3} {r["logL"]:>12,.2f} {r["logL"] / n_bins:>+9.4f} '
                  f'{delta:>11} {pk:>4} {aic:>12,.1f} {bic:>12,.1f} '
                  f'{r["n_iter_run"]:>6}{flag}')

        # --- monotonicity ---
        lls = [r['logL'] for r in rows]
        drops = [(rows[i]['k'], rows[i + 1]['k'], lls[i + 1] - lls[i])
                 for i in range(len(lls) - 1) if lls[i + 1] < lls[i]]
        if drops:
            print('  [!] logL DECREASED with k (impossible at the global optimum):')
            for a, b, d in drops:
                print(f'      k={a} -> k={b}: {d:+,.2f}')
        else:
            print('  logL monotone non-decreasing in k: OK')

        # --- sign check ---
        if true_ll is not None:
            exp_neg = meta.get('expect_negative_aic')
            got_neg = (2 * n_params(best_aic['k']) - 2 * best_aic['logL']) < 0
            verdict = 'as predicted' if exp_neg == got_neg else 'MISMATCH'
            print(f'  AIC sign: {"negative" if got_neg else "positive"} '
                  f'(predicted {"negative" if exp_neg else "positive"}) -- {verdict}')

        # --- parameter recovery at k_true ---
        if k_true:
            row = next((r for r in rows if r['k'] == k_true), None)
            if row and meta.get('means'):
                tm = np.asarray(meta['means'])
                ts = np.asarray(meta['sigmas']) ** 2
                fm = np.asarray(sorted(row['means']))
                fv = np.asarray(row['variances'])
                print(f'  recovery at k={k_true}:')
                print(f'    means  true {np.round(tm, 4).tolist()}  '
                      f'fitted {np.round(fm, 4).tolist()}')
                print(f'    vars   true {np.round(ts, 4).tolist()}  '
                      f'fitted {np.round(fv, 4).tolist()}')
                if len(fm) == len(tm):
                    print(f'    max |mean error| = {np.abs(fm - tm).max():.4f}')


if __name__ == '__main__':
    main()
