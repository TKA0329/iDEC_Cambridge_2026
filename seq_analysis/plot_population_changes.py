#!/usr/bin/env python3
"""Compare count-weighted sequence-property populations across 2-3 samples.

Each input CSV must contain a sequence column, a nonnegative integer count
column, and numeric property columns. The script writes count-weighted means,
medians, population standard deviations, and pairwise Mann–Whitney U tests,
then plots the means as bars and marks significant comparisons.
"""

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import norm

matplotlib.use('Agg')
import matplotlib.pyplot as plt


SEQUENCE_CANDIDATES = ('sequence', 'aa_variant')
NON_PROPERTY_COLUMNS = {
    'id', 'variant_id', 'motif1_sequence', 'error', 'stable',
    'helix_breaker_flag', 'cluster', 'fcr', 'ncpr'
}
COLORS = ['#9fc8c8', '#65a6a6', '#176f6f']
PROPERTY_TITLES = {
    'pred_cv_JperK': 'Heat capacity (J/K)',
    'pred_rog_A': 'Radius of gyration (Å)',
    'pred_tau_fs': 'End-to-end decorrelation time (fs)',
}


def resolve_column(frame, requested, candidates, role, path):
    """Resolve a requested or conventional column name case-insensitively."""
    lookup = {str(column).lower(): column for column in frame.columns}
    if requested:
        found = lookup.get(requested.lower())
        if found is None:
            raise ValueError(f'{path}: missing {role} column {requested!r}')
        return found
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    raise ValueError(f'{path}: could not detect a {role} column')


def load_input(path, sequence_col, count_col):
    """Load one input and validate sequence and count values."""
    frame = pd.read_csv(path)
    seq = resolve_column(frame, sequence_col, SEQUENCE_CANDIDATES, 'sequence', path)
    count = resolve_column(frame, count_col, ('count',), 'count', path)

    if frame[seq].isna().any():
        raise ValueError(f'{path}: missing sequences')
    frame[seq] = frame[seq].astype(str).str.strip()
    if frame[seq].eq('').any():
        raise ValueError(f'{path}: empty sequences')

    frame[count] = pd.to_numeric(frame[count], errors='raise')
    values = frame[count].to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or (values % 1 != 0).any():
        raise ValueError(f'{path}: counts must be finite nonnegative integers')
    frame[count] = frame[count].astype('int64')
    if frame[count].sum() == 0:
        raise ValueError(f'{path}: total count must be positive')
    return frame, seq, count


def numeric_property(frame, column):
    """Return numeric values, or None when a column is not a usable property."""
    if pd.api.types.is_bool_dtype(frame[column]):
        return None
    try:
        values = pd.to_numeric(frame[column], errors='raise')
    except (TypeError, ValueError):
        return None
    array = values.to_numpy(dtype=float)
    if np.isinf(array).any() or not np.isfinite(array).any():
        return None
    return values


def select_features(loaded, requested):
    """Select explicit features or auto-detect shared numeric properties."""
    frames = [item[0] for item in loaded]
    excluded = set(NON_PROPERTY_COLUMNS)
    for _, seq, count in loaded:
        excluded.update((str(seq).lower(), str(count).lower()))

    if requested:
        features = [value.strip() for value in requested.split(',')]
        if not features or not all(features) or len(features) != len(set(features)):
            raise ValueError('--features must contain unique nonempty names')
        for feature in features:
            for frame in frames:
                if feature not in frame:
                    raise ValueError(f'Missing requested property {feature!r}')
                if numeric_property(frame, feature) is None:
                    raise ValueError(f'Property {feature!r} must be numeric and non-boolean')
        return features

    shared = set(frames[0].columns)
    for frame in frames[1:]:
        shared.intersection_update(frame.columns)
    features = []
    for column in frames[0].columns:
        if column not in shared or str(column).lower() in excluded:
            continue
        if all(numeric_property(frame, column) is not None for frame in frames):
            features.append(column)
    if not features:
        raise ValueError('No shared numeric property columns were detected')
    return features


def expanded_population_median(values, weights):
    """Median obtained if each value were repeated by its integer count."""
    order = np.argsort(values, kind='stable')
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    total = int(cumulative[-1])
    lower_position = (total - 1) // 2
    upper_position = total // 2
    lower = values[np.searchsorted(cumulative, lower_position, side='right')]
    upper = values[np.searchsorted(cumulative, upper_position, side='right')]
    return float((lower + upper) / 2)


def summarise(loaded, labels, paths, features):
    """Calculate count-weighted population statistics for every property."""
    rows = []
    for (frame, seq_col, count_col), label, path in zip(loaded, labels, paths):
        counts = frame[count_col].to_numpy(dtype=np.int64)
        total_reads = int(counts.sum())
        for feature in features:
            values = pd.to_numeric(frame[feature], errors='raise').to_numpy(dtype=float)
            valid = np.isfinite(values) & (counts > 0)
            valid_values = values[valid]
            valid_counts = counts[valid]
            reads_used = int(valid_counts.sum())
            if reads_used == 0:
                raise ValueError(f'{path}: {feature} has no finite values with positive counts')
            mean = float(np.average(valid_values, weights=valid_counts))
            population_sd = float(
                np.sqrt(np.average((valid_values - mean) ** 2, weights=valid_counts))
            )
            rows.append(
                {
                    'sample': label,
                    'input_file': str(path),
                    'property': feature,
                    'total_reads': total_reads,
                    'reads_used': reads_used,
                    'reads_excluded_missing_property': total_reads - reads_used,
                    'sequence_rows_used': int(valid.sum()),
                    'unique_sequences_used': int(frame.loc[valid, seq_col].nunique()),
                    'population_mean': mean,
                    'population_median': expanded_population_median(
                        valid_values, valid_counts
                    ),
                    'population_sd': population_sd,
                }
            )
    return pd.DataFrame(rows)


def weighted_mann_whitney(values_a, weights_a, values_b, weights_b):
    """Two-sided asymptotic Mann–Whitney U test for frequency-weighted data."""
    valid_a = np.isfinite(values_a) & (weights_a > 0)
    valid_b = np.isfinite(values_b) & (weights_b > 0)
    values_a, weights_a = values_a[valid_a], weights_a[valid_a]
    values_b, weights_b = values_b[valid_b], weights_b[valid_b]
    n_a, n_b = int(weights_a.sum()), int(weights_b.sum())
    if n_a == 0 or n_b == 0:
        raise ValueError('Mann–Whitney U requires positive counts in both samples')

    unique_b, inverse_b = np.unique(values_b, return_inverse=True)
    grouped_b = np.bincount(inverse_b, weights=weights_b)
    cumulative_b = np.concatenate(([0.0], np.cumsum(grouped_b)))
    left = np.searchsorted(unique_b, values_a, side='left')
    right = np.searchsorted(unique_b, values_a, side='right')
    less = cumulative_b[left]
    equal = cumulative_b[right] - cumulative_b[left]
    u_a = float(np.dot(weights_a, less + 0.5 * equal))

    pooled_values = np.concatenate((values_a, values_b))
    pooled_weights = np.concatenate((weights_a, weights_b))
    _, pooled_inverse = np.unique(pooled_values, return_inverse=True)
    ties = np.bincount(pooled_inverse, weights=pooled_weights)
    total = n_a + n_b
    mean_u = n_a * n_b / 2
    tie_term = np.sum(ties ** 3 - ties)
    variance_u = n_a * n_b / 12 * (
        total + 1 - tie_term / (total * (total - 1))
    )
    if variance_u <= 0:
        p_value = 1.0
    else:
        # Continuity correction matches the usual asymptotic two-sided test.
        distance = max(0.0, abs(u_a - mean_u) - 0.5)
        p_value = float(2 * norm.sf(distance / np.sqrt(variance_u)))
    return u_a, p_value, n_a, n_b


def holm_adjust(p_values):
    """Holm family-wise error correction, returned in original order."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * p_values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def mann_whitney_tests(loaded, labels, features, alpha, adjustment):
    """Run all pairwise tests for every plotted property."""
    rows = []
    for feature in features:
        for first_index in range(len(loaded) - 1):
            frame_a, _, count_a = loaded[first_index]
            for second_index in range(first_index + 1, len(loaded)):
                frame_b, _, count_b = loaded[second_index]
                u_statistic, p_value, n_a, n_b = weighted_mann_whitney(
                    pd.to_numeric(frame_a[feature], errors='raise').to_numpy(dtype=float),
                    frame_a[count_a].to_numpy(dtype=np.int64),
                    pd.to_numeric(frame_b[feature], errors='raise').to_numpy(dtype=float),
                    frame_b[count_b].to_numpy(dtype=np.int64),
                )
                rows.append(
                    {
                        'property': feature,
                        'sample_a': labels[first_index],
                        'sample_b': labels[second_index],
                        'population_size_a': n_a,
                        'population_size_b': n_b,
                        'u_statistic_a': u_statistic,
                        'p_value': p_value,
                    }
                )
    tests = pd.DataFrame(rows)
    tests['adjusted_p_value'] = (
        holm_adjust(tests.p_value) if adjustment == 'holm' else tests.p_value
    )
    tests['significant'] = tests.adjusted_p_value < alpha
    return tests


def property_title(feature):
    """Return a publication-friendly property title."""
    return PROPERTY_TITLES.get(feature, str(feature).replace('_', ' '))


def property_filename(feature):
    """Return a filesystem-safe, recognisable filename for one property."""
    safe = re.sub(r'[^A-Za-z0-9._-]+', '_', str(feature)).strip('._')
    return f'{safe or "property"}.png'


def format_p_value(value):
    """Compact p-value label for plot annotations."""
    return f'{value:.2e}' if value < 0.001 else f'{value:.3f}'


def draw_property(ax, data, tests, feature, labels, whiskers):
    """Draw one property's count-weighted means and optional SD whiskers."""
    x = np.arange(len(labels))
    data = data.set_index('sample').loc[labels]
    means = data.population_mean.to_numpy()
    errors = data.population_sd.to_numpy() if whiskers == 'sd' else None
    ax.bar(
        x,
        means,
        yerr=errors,
        capsize=5 if errors is not None else 0,
        color=COLORS[:len(labels)],
        edgecolor='#333333',
        linewidth=0.6,
    )
    ax.set_xticks(
        x,
        labels,
        rotation=15 if any(len(label) > 12 for label in labels) else 0,
    )
    ax.set_title(property_title(feature))
    ax.set_ylabel('Count-weighted population mean')
    ax.grid(axis='y', alpha=0.2)
    ax.spines[['top', 'right']].set_visible(False)
    significant = tests[(tests.property == feature) & tests.significant]
    if not significant.empty:
        lines = [
            f'* {row.sample_a} vs {row.sample_b}: p_adj={format_p_value(row.adjusted_p_value)}'
            for row in significant.itertuples()
        ]
        ax.text(
            0.02,
            0.98,
            '\n'.join(lines),
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=8,
            fontweight='bold',
            bbox={'facecolor': 'white', 'alpha': 0.75, 'edgecolor': 'none'},
        )


def plot_summary(summary, tests, features, labels, whiskers, title, outdir):
    """Draw combined and individual mean bar charts for every property."""
    ncols = min(3, len(features))
    nrows = math.ceil(len(features) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.6 * ncols, 3.6 * nrows + 0.7),
        squeeze=False,
    )
    for ax, feature in zip(axes.flat, features):
        draw_property(
            ax,
            summary[summary.property == feature],
            tests,
            feature,
            labels,
            whiskers,
        )
    for ax in axes.flat[len(features):]:
        ax.remove()

    note = 'Bars: count-weighted means.'
    if whiskers == 'sd':
        note += ' Whiskers: ±1 count-weighted population SD (not confidence intervals).'
    note += ' * indicates a significant pairwise Mann–Whitney U test after correction.'
    fig.suptitle(title, fontsize=14)
    fig.text(0.5, 0.01, note, ha='center', fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 0.96))
    for extension in ('png', 'pdf', 'svg'):
        fig.savefig(outdir / f'population_property_means.{extension}', dpi=300)
    plt.close(fig)

    individual_dir = outdir / 'individual_plots'
    individual_dir.mkdir(parents=True, exist_ok=True)
    for feature in features:
        individual_fig, individual_ax = plt.subplots(figsize=(5.4, 4.2))
        draw_property(
            individual_ax,
            summary[summary.property == feature],
            tests,
            feature,
            labels,
            whiskers,
        )
        individual_fig.tight_layout()
        individual_fig.savefig(
            individual_dir / property_filename(feature),
            dpi=300,
        )
        plt.close(individual_fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, nargs='+')
    parser.add_argument('--labels', nargs='+', help='One display label per input')
    parser.add_argument('--features', help='Comma-separated property names; default: auto-detect')
    parser.add_argument('--sequence-col', help='Sequence column; default: sequence or aa_variant')
    parser.add_argument('--count-col', help='Count column; default: count')
    parser.add_argument('--whiskers', choices=('sd', 'none'), default='sd')
    parser.add_argument(
        '--alpha', type=float, default=0.05,
        help='Significance threshold for adjusted p-values (default: 0.05)',
    )
    parser.add_argument(
        '--p-adjust', choices=('holm', 'none'), default='holm',
        help='Multiple-testing correction across all properties and pairs',
    )
    parser.add_argument('--title', default='Population property changes')
    parser.add_argument('--outdir', type=Path, default=Path('population_changes'))
    args = parser.parse_args()

    if len(args.inputs) not in (2, 3):
        parser.error('Provide exactly two or three input CSVs')
    if not 0 < args.alpha < 1:
        parser.error('--alpha must be between 0 and 1')
    default_labels = ['Before', 'After'] if len(args.inputs) == 2 else [
        'Before', 'Intermediate', 'After'
    ]
    labels = args.labels or default_labels
    if len(labels) != len(args.inputs) or len(set(labels)) != len(labels):
        parser.error('--labels must contain one unique label per input')

    loaded = [
        load_input(path, args.sequence_col, args.count_col)
        for path in args.inputs
    ]
    features = select_features(loaded, args.features)
    summary = summarise(loaded, labels, args.inputs, features)
    tests = mann_whitney_tests(
        loaded, labels, features, args.alpha, args.p_adjust
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.outdir / 'population_statistics.csv', index=False)
    tests.to_csv(args.outdir / 'mann_whitney_tests.csv', index=False)
    plot_summary(
        summary, tests, features, labels, args.whiskers, args.title, args.outdir
    )
    metadata = {
        'input_files': [str(path) for path in args.inputs],
        'labels': labels,
        'properties': features,
        'whiskers': args.whiskers,
        'mann_whitney_alpha': args.alpha,
        'mann_whitney_p_adjustment': args.p_adjust,
        'multiple_testing_family': 'all properties and pairwise sample comparisons',
        'mann_whitney_interpretation': 'rank/distribution difference, not specifically a mean test',
        'weighting': 'sequence row weighted by count',
        'standard_deviation': 'population (ddof=0)',
        'median': 'median of the count-expanded population',
    }
    (args.outdir / 'run_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Analysed {len(features)} properties across {len(args.inputs)} samples')
    print(
        f'{int(tests.significant.sum())} of {len(tests)} pairwise tests were '
        f'significant at adjusted p < {args.alpha:g}'
    )
    print(f'Saved statistics, tests, combined plots, and individual PNGs to {args.outdir}')



if __name__ == '__main__':
    main()
