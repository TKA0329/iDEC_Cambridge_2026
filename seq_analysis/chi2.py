"""
Monte Carlo simulation of chi-square test for comparing populations of variants before and after selection.

To use:

Input columns: aa_variant, count (overrideable with options --variant-col and --count-col).

"""

import argparse
import csv
import hashlib
import json

from collections import Counter
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import beta, random_table


# returns a counter (dict-like object) of {variant: count} from a CSV file
def read_counts(path, variant_col, count_col):
    counts = Counter()
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not {variant_col, count_col}.issubset(reader.fieldnames):
            raise ValueError(f'{path}: required columns: {variant_col}, {count_col}')
        for line, row in enumerate(reader, 2):
            seq = (row.get(variant_col) or '').strip().upper()
            raw = (row.get(count_col) or '').strip()
            if not seq or any(c.isspace() for c in seq):
                raise ValueError(f'{path}:{line}: blank variant or internal whitespace')
            # Reject fractional, negative and missing counts; do not silently round.
            if not raw.isascii() or not raw.isdigit():
                raise ValueError(f'{path}:{line}: count must be a nonnegative integer')
            counts[seq] += int(raw)
    if sum(counts.values()) == 0:
        raise ValueError(f'{path}: Sample has no positive counts')
    return counts

def monte_carlo_test(observed, simulations=100000, seed=42, batch_size=1000):
    """
    Null hypothesis: the selected vs unselected populations are the same
    If null hypothesis is true, the observed contingency table is one of many 
    possible tables with the same row and column totals. The Monte Carlo simulation generates 
    random tables with the same margins and computes the chi-square statistic for each.
    The p-value is estimated as the proportion of simulated statistics that are greater than
    or equal to the observed statistic.

    Return statistic, expected counts, simulated-tail count and corrected p.

    random_table samples from the conditional independence distribution;
    it does NOT give every possible contingency table equal probability.

    Returns:
        statistic: Pearson chi-square statistic for the observed table
        expected: expected counts under the null hypothesis
        exceedances: number of simulated statistics >= observed statistic
        p_value: Monte Carlo p-value with plus-one correction
        [lo, hi]: 95% Clopper-Pearson interval
    """
    if simulations < 1 or batch_size < 1:
        raise ValueError('Simulations and batch size must be positive')
    row_totals = observed.sum(axis=1)
    col_totals = observed.sum(axis=0)
    expected = np.outer(row_totals.astype(float), col_totals) / observed.sum()

    # 1. Pearson chi-square statistic, without Yates correction.
    statistic = float(np.sum((observed - expected) ** 2 / expected))
    rng = np.random.default_rng(seed)
    exceedances = 0 # track how many simulations exceed the observed statistic
    # Include numerical ties in the upper tail.
    tolerance = 100 * np.finfo(float).eps * max(1.0, abs(statistic))
    for start in range(0, simulations, batch_size):
        size = min(batch_size, simulations - start)
        # 2. Null tables preserve BOTH sample sizes and pooled variant totals.
        tables = random_table.rvs(row_totals, col_totals, size=size,
                                  method='patefield', random_state=rng)
        # 3. Same expected counts for every table because margins are fixed.
        simulated_chi2 = np.sum((tables - expected) ** 2 / expected, axis=(1, 2))
        exceedances += int(np.count_nonzero(simulated_chi2 >= statistic - tolerance))
    # 4. Plus-one Monte Carlo p-value: never report zero.
    p_value = (exceedances + 1) / (simulations + 1)
    # Binomial interval describes Monte Carlo numerical uncertainty ONLY.
    lo = 0.0 if exceedances == 0 else float(beta.ppf(
        0.025, exceedances, simulations - exceedances + 1))
    hi = 1.0 if exceedances == simulations else float(beta.ppf(
        0.975, exceedances + 1, simulations - exceedances))
    return statistic, expected, exceedances, p_value, [lo, hi]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('unselected', type=Path)
    parser.add_argument('selected', type=Path)
    parser.add_argument('--variant-col', default='aa_variant', help='Column name for the variant identifier')
    parser.add_argument('--count-col', default='count', help='Column name for the read count')
    parser.add_argument('--simulations', type=int, default=100000, help='Number of Monte Carlo simulations to run. More simulations give more accurate p-values but take longer.')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    parser.add_argument('--batch-size', type=int, default=1000, help='Number of simulations to run in each batch. Larger batches are faster but use more memory.')
    parser.add_argument('--outdir', type=Path, default=Path('chi2_results'), help='Output directory for results')
    args = parser.parse_args()
    try:
        baseline = read_counts(args.unselected, args.variant_col, args.count_col)
        selected = read_counts(args.selected, args.variant_col, args.count_col)
        variants = sorted(v for v in baseline.keys() | selected.keys()
                          if baseline[v] + selected[v] > 0)
        if len(variants) < 2:
            raise ValueError('At least two positive pooled variant categories are required')
        total = sum(baseline.values()) + sum(selected.values())
        if total > np.iinfo(np.int64).max:
            raise ValueError('Read total exceeds supported integer range')
        observed = np.array([[baseline[v] for v in variants],
                             [selected[v] for v in variants]], dtype=np.int64)
        stat, expected, hits, p, interval = monte_carlo_test(
            observed, args.simulations, args.seed, args.batch_size)
    except ValueError as exc:
        parser.error(str(exc))
    args.outdir.mkdir(parents=True, exist_ok=True)
    contributions = np.sum((observed - expected) ** 2 / expected, axis=0)
    with open(args.outdir / 'aligned_counts.csv', 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['aa_variant', 'unselected_count', 'selected_count',
                         'unselected_frequency', 'selected_frequency',
                         'unselected_expected', 'selected_expected', 'chi2_contribution'])
        for i, variant in enumerate(variants):
            writer.writerow([variant, int(observed[0, i]), int(observed[1, i]),
                             observed[0, i]/observed[0].sum(),
                             observed[1, i]/observed[1].sum(),
                             expected[0, i], expected[1, i], contributions[i]])
    summary = {
        'test': 'Pearson chi-square with fixed-margin Monte Carlo calibration',
        'rows': ['unselected', 'selected'],
        'input_files': [str(args.unselected), str(args.selected)],
        'input_sha256': [hashlib.sha256(pth.read_bytes()).hexdigest()
                         for pth in (args.unselected, args.selected)],
        'sample_read_totals': observed.sum(axis=1).tolist(),
        'observed_variants_per_sample': (observed > 0).sum(axis=1).tolist(),
        'pooled_variants': len(variants),
        'pearson_chi2': stat,
        'nominal_degrees_of_freedom': len(variants)-1,
        'expected_cells_below_5': int((expected < 5).sum()),
        'total_cells': int(expected.size),
        'minimum_expected_count': float(expected.min()),
        'simulations': args.simulations, 'seed': args.seed,
        'batch_size': args.batch_size, 'sampler': 'patefield',
        'simulated_statistics_at_least_observed': hits,
        'monte_carlo_p_plus_one': p,
        'monte_carlo_tail_probability_95pct_binomial_interval': interval,
        'cramers_v_descriptive': float(np.sqrt(stat/observed.sum())),
        'numpy_version': np.__version__, 'scipy_version': scipy.__version__,
        'interpretation': ('Read-composition comparison conditional on acceptance rules and '
                           'independent representative read sampling. Not a test of '
                           'reproducibility across biological selection replicates.'),
    }
    (args.outdir / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))
    if hits == 0:
        print('No simulated statistic reached the observed value. The p-value is at '
              'the simulation resolution floor, not an accurately resolved extreme tail.')


if __name__ == '__main__':
    main()
