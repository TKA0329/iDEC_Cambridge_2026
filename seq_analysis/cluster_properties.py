#!/usr/bin/env python3
"""Compare amino-acid sequence properties across two or three samples.

The script pools the unique sequences from all samples, clusters them using
standardised biochemical properties, and projects the same property data onto
principal components for plotting. Read counts affect point sizes, cluster
abundances, and weighted centres; they do *not* affect the fitted clusters.

Typical usage::

    python cluster_properties.py untreated.csv intermediate.csv selected.csv --outdir results

Dependencies: numpy, pandas, scipy, scikit-learn, and matplotlib.
"""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scipy, sklearn
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, cut_tree, dendrogram
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

# Use a non-interactive backend so the script works on servers and in pipelines.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Default numerical properties used to describe each amino-acid sequence.
DEFAULT_FEATURES = 'gravy,charge_at_pH7,instability_index,turn_fraction,helix_sheet_sum,pred_cv_JperK,pred_rog_A,pred_tau_fs'


def load(path, features, seq_col, count_col):
    """Load one sample CSV and validate the columns needed by the analysis.

    Rows with a count of zero are removed because they contribute no reads to
    either the abundance calculations or the plots. Duplicate sequence rows
    are retained here and summed later, after both samples have been loaded.
    """
    d = pd.read_csv(path)
    required = [seq_col, count_col] + features
    missing = set(required) - set(d.columns)

    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}')

    # Upstream analyser failures make the property columns unreliable.
    if 'error' in d and d.error.fillna('').astype(str).str.strip().ne('').any():
        raise ValueError(f'{path}: analyser errors present; resolve before clustering')

    if d[seq_col].isna().any():
        raise ValueError('Missing sequence')

    # Normalise sequence text so case and surrounding whitespace do not create
    # false "different" variants.
    d[seq_col] = d[seq_col].astype(str).str.strip().str.upper()
    if d[seq_col].eq('').any():
        raise ValueError('Empty sequence')

    # Fail early on non-numeric, infinite, or otherwise unusable values.
    for f in [count_col] + features:
        d[f] = pd.to_numeric(d[f], errors='raise')
        if not np.isfinite(d[f]).all():
            raise ValueError(f'Nonfinite values in {f}')

    if (d[count_col] < 0).any() or (d[count_col] % 1 != 0).any():
        raise ValueError('Counts must be nonnegative integers')

    d[count_col] = d[count_col].astype('int64')
    return d.loc[d[count_col] > 0, required]


def main():
    """Parse command-line options and run the complete comparison workflow."""

    # ------------------------------------------------------------------
    # 1. Parse and validate command-line arguments.
    # ------------------------------------------------------------------
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        'inputs',
        type=Path,
        nargs='+',
        help='Two or three sample CSVs in experimental order',
    )
    p.add_argument(
        '--sample-labels',
        nargs='+',
        help='Optional display labels in the same order as the input CSVs',
    )
    p.add_argument(
        '--features',
        default=DEFAULT_FEATURES,
        help='Comma-separated numeric features; defaults tailored to SSB/NEXT',
    )
    p.add_argument(
        '--clusters',
        default='auto',
        help='auto or integer 2-9; auto maximises silhouette over 2-9',
    )
    p.add_argument('--sequence-col', default='aa_variant')
    p.add_argument('--count-col', default='count')
    p.add_argument(
        '--min-pooled-count',
        type=int,
        default=1,
        help='Use 2 for separate singleton-exclusion sensitivity analysis',
    )
    p.add_argument('--outdir', type=Path, default=Path('cluster_results'))

    a = p.parse_args()
    fs = [feature.strip() for feature in a.features.split(',')]
    if len(fs) != len(set(fs)) or not all(fs):
        p.error('Features must be nonempty and unique')
    if a.min_pooled_count < 1:
        p.error('Minimum pooled count must be >=1')
    if len(a.inputs) not in (2, 3):
        p.error('Provide exactly two or three input CSVs')

    sample_names = (
        ['unselected', 'selected']
        if len(a.inputs) == 2
        else ['unselected', 'intermediate', 'selected']
    )
    default_labels = (
        ['Untreated', 'Selected']
        if len(a.inputs) == 2
        else ['Untreated', 'Intermediate', 'Selected']
    )
    sample_labels = a.sample_labels or default_labels
    if len(sample_labels) != len(a.inputs):
        p.error('--sample-labels must provide one label per input CSV')

    # ------------------------------------------------------------------
    # 2. Load all samples and align them by amino-acid sequence.
    # ------------------------------------------------------------------
    sample_data = {
        name: load(path, fs, a.sequence_col, a.count_col)
        for name, path in zip(sample_names, a.inputs)
    }
    combined = pd.concat(sample_data.values(), ignore_index=True)

    # A sequence's calculated traits should be intrinsic to that sequence.
    # Check that duplicate rows agree, both within and across samples.
    for seq, group in combined.groupby(a.sequence_col):
        values = group[fs].to_numpy()
        first_values = group[fs].iloc[0].to_numpy()
        if not np.allclose(values, first_values, rtol=1e-8, atol=1e-8):
            raise ValueError(f'Conflicting feature values for sequence {seq}')

    # `traits` has one row per unique sequence. `counts` contains the number of
    # reads for that sequence in each sample; an absent sequence receives zero.
    traits = combined.groupby(a.sequence_col, sort=True)[fs].first()
    counts = pd.DataFrame(
        {
            name: data.groupby(a.sequence_col)[a.count_col].sum()
            for name, data in sample_data.items()
        }
    ).fillna(0).astype(int)
    counts = counts.reindex(traits.index)

    # Optionally remove sequences whose total reads across all samples are
    # below the threshold. The default of 1 retains every observed sequence.
    original_totals = counts.sum()
    keep = counts.sum(axis=1) >= a.min_pooled_count
    excluded = counts.loc[~keep].copy()
    counts = counts.loc[keep]
    traits = traits.loc[keep]

    if len(traits) < 3 or (counts.sum() == 0).any():
        raise ValueError('Need >=3 sequences and positive totals in every sample')

    # Constant properties contain no clustering information and StandardScaler
    # cannot give them meaningful variance, so exclude them from the model.
    active = [f for f in fs if traits[f].std(ddof=0) > 1e-12]
    dropped = [f for f in fs if f not in active]
    if not active:
        raise ValueError('All selected features are constant')

    # Standardisation gives every property mean 0 and variance 1. Without it,
    # properties measured on larger numerical scales would dominate distances.
    scaler = StandardScaler()
    X = scaler.fit_transform(traits[active])

    # ------------------------------------------------------------------
    # 3. Fit Ward hierarchical clustering in the full property space.
    # ------------------------------------------------------------------
    # There is one point per unique sequence. Counts and sample labels do not
    # influence the cluster fit. Ward linkage successively merges the pair of
    # groups causing the smallest increase in within-cluster variance.
    Z = linkage(X, method='ward', metric='euclidean', optimal_ordering=True)

    # Score each allowed cluster count. A higher silhouette score means that
    # sequences are, on average, closer to their own cluster than other ones.
    diagnostic = []
    for candidate_k in range(2, min(9, len(X) - 1) + 1):
        candidate_labels = cut_tree(Z, n_clusters=candidate_k).ravel()
        diagnostic.append(
            {
                'k': candidate_k,
                'silhouette': silhouette_score(X, candidate_labels),
                'smallest_cluster': int(np.bincount(candidate_labels).min()),
            }
        )
    diag = pd.DataFrame(diagnostic)

    if a.clusters == 'auto':
        # Break equal-score ties in favour of the smaller, simpler solution.
        best = diag.sort_values(['silhouette', 'k'], ascending=[False, True]).iloc[0]
        k = int(best.k)
    else:
        k = int(a.clusters)
        if not 2 <= k <= min(9, len(X) - 1):
            raise ValueError(
                'Requested cluster count must be 2-9 and less than number of sequences'
            )

    # scipy numbers clusters from zero; add one for user-facing labels C1, C2, ...
    labels = cut_tree(Z, n_clusters=k).ravel() + 1

    # ------------------------------------------------------------------
    # 4. Save clustering diagnostics and calculate PCA coordinates.
    # ------------------------------------------------------------------
    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    diag.to_csv(out / 'cluster_diagnostics.csv', index=False)
    traits.corr().to_csv(out / 'feature_correlations.csv')
    pd.DataFrame(
        {
            'feature': active,
            'pooled_mean': scaler.mean_,
            'pooled_sd': scaler.scale_,
        }
    ).to_csv(out / 'scaling.csv', index=False)
    pd.DataFrame(
        Z,
        columns=['left', 'right', 'ward_distance', 'members'],
    ).to_csv(out / 'linkage.csv', index=False)
    excluded.to_csv(out / 'excluded_counts.csv')

    # PCA rotates the standardised feature space onto axes of decreasing
    # variance. It is ONLY a 2-D display projection; clustering used all active
    # properties in X above. A third component is retained in the CSV if valid.
    ncomp = min(3, len(active), len(X))
    pca = PCA(n_components=ncomp, svd_solver='full')
    coords = pca.fit_transform(X)
    if ncomp < 2:
        # Matplotlib still needs a y-coordinate when only one property varies.
        coords = np.column_stack([coords, np.zeros(len(coords))])

    # Create the per-sequence results table: cluster, counts, frequencies, raw
    # properties, and PCA coordinates all remain linked by the sequence index.
    members = traits.copy()
    members.insert(0, 'variant_id', [f'V{i + 1:03d}' for i in range(len(traits))])
    members.insert(1, 'cluster', labels)
    for s in counts:
        members[s + '_reads'] = counts[s]
        members[s + '_frequency'] = counts[s] / counts[s].sum()
    for i in range(coords.shape[1]):
        members[f'PC{i + 1}'] = coords[:, i]
    members.to_csv(out / 'members.csv')
    pd.DataFrame(
        pca.components_.T,
        index=active,
        columns=[f'PC{i + 1}' for i in range(ncomp)],
    ).to_csv(out / 'pca_loadings.csv')
    # ------------------------------------------------------------------
    # 5. Summarise abundance and properties within each cluster.
    # ------------------------------------------------------------------
    summary = []
    trait_rows = []
    for c in range(1, k + 1):
        mask = labels == c
        cluster_traits = traits.loc[mask]
        cluster_counts = counts.loc[mask]

        # Count both distinct sequences and total reads. The fractions show the
        # share of all reads assigned to this cluster in each sample.
        row = {'cluster': c, 'unique_sequences_pooled': int(mask.sum())}
        for s in counts:
            row[s + '_unique_sequences'] = int((cluster_counts[s] > 0).sum())
            row[s + '_reads'] = int(cluster_counts[s].sum())
            row[s + '_fraction'] = float(
                cluster_counts[s].sum() / counts[s].sum()
            )
        # Retain the original overall change column and add adjacent changes
        # when an intermediate sample is present.
        row['percentage_point_change'] = 100 * (
            row['selected_fraction'] - row['unselected_fraction']
        )
        for previous, current in zip(sample_names, sample_names[1:]):
            row[f'{current}_minus_{previous}_percentage_points'] = 100 * (
                row[f'{current}_fraction'] - row[f'{previous}_fraction']
            )
        summary.append(row)

        # Report each property's ordinary distribution across unique sequences,
        # plus a read-weighted distribution for each individual sample.
        for f in fs:
            property_summary = {
                'cluster': c,
                'feature': f,
                'unique_sequence_sum': float(cluster_traits[f].sum()),
                'unique_sequence_mean': float(cluster_traits[f].mean()),
                'unique_sequence_median': float(cluster_traits[f].median()),
                'unique_sequence_sd': float(cluster_traits[f].std(ddof=0)),
                'min': float(cluster_traits[f].min()),
                'max': float(cluster_traits[f].max()),
            }
            for s in counts:
                weights = cluster_counts[s].to_numpy()
                total = weights.sum()
                weighted_mean = (
                    float(np.average(cluster_traits[f], weights=weights))
                    if total
                    else np.nan
                )
                property_summary[s + '_read_weighted_sum'] = float(
                    np.dot(weights, cluster_traits[f])
                )
                property_summary[s + '_read_weighted_mean'] = weighted_mean
                property_summary[s + '_read_weighted_sd'] = (
                    float(
                        np.sqrt(
                            np.average(
                                (cluster_traits[f] - weighted_mean) ** 2,
                                weights=weights,
                            )
                        )
                    )
                    if total
                    else np.nan
                )
            trait_rows.append(property_summary)

    summary = pd.DataFrame(summary)
    summary.to_csv(out / 'cluster_counts.csv', index=False)
    pd.DataFrame(trait_rows).to_csv(out / 'cluster_traits.csv', index=False)

    # These tables put cluster read counts into formats convenient for later
    # contingency-table or per-cluster pairwise Fisher exact tests.
    read_columns = [f'{sample}_reads' for sample in sample_names]
    contingency = summary.set_index('cluster')[read_columns].T
    contingency.to_csv(out / 'cluster_contingency.csv')

    if len(sample_names) == 2:
        # Preserve the original two-sample export schema.
        fisher_rows = [
            {
                'cluster': int(row.cluster),
                'selected_in': int(row.selected_reads),
                'selected_out': int(counts.selected.sum() - row.selected_reads),
                'unselected_in': int(row.unselected_reads),
                'unselected_out': int(counts.unselected.sum() - row.unselected_reads),
            }
            for row in summary.itertuples()
        ]
    else:
        # Three samples produce one 2x2 input table for every sample pair.
        fisher_rows = []
        for row in summary.itertuples():
            for first_index, first in enumerate(sample_names[:-1]):
                for second in sample_names[first_index + 1:]:
                    first_in = int(getattr(row, f'{first}_reads'))
                    second_in = int(getattr(row, f'{second}_reads'))
                    fisher_rows.append(
                        {
                            'cluster': int(row.cluster),
                            'sample_a': first,
                            'sample_b': second,
                            'sample_a_in': first_in,
                            'sample_a_out': int(counts[first].sum() - first_in),
                            'sample_b_in': second_in,
                            'sample_b_out': int(counts[second].sum() - second_in),
                        }
                    )
    pd.DataFrame(fisher_rows).to_csv(out / 'fisher_input_tables.csv', index=False)

    # ------------------------------------------------------------------
    # 6. Plot the sample comparison and clustering diagnostics.
    # ------------------------------------------------------------------
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    colors = plt.get_cmap('tab10')
    pcvar = pca.explained_variance_ratio_

    # All panels use identical axes and point-size scaling. A point represents
    # a unique sequence; its area represents its fraction of reads in that
    # sample. Colour represents its pooled Ward cluster.
    fig, axes = plt.subplots(
        1,
        len(sample_names),
        figsize=(5.3 * len(sample_names), 4.8),
        sharex=True,
        sharey=True,
    )
    for ax, sample, title in zip(axes, sample_names, sample_labels):
        sample_counts = counts[sample].to_numpy()
        for c in range(1, k + 1):
            visible = (labels == c) & (sample_counts > 0)
            point_areas = sample_counts[visible] / counts[sample].sum() * 1500
            ax.scatter(
                coords[visible, 0],
                coords[visible, 1],
                s=point_areas,
                color=colors(c - 1),
                alpha=0.75,
                edgecolors='black',
                linewidth=0.4,
                label=f'C{c}',
            )

            # The cross is the cluster's centre in PCA space, weighted by read
            # counts in this sample. It can move between samples even though
            # cluster membership itself is fixed across all panels.
            weights = sample_counts[visible]
            if weights.sum():
                center = np.average(coords[visible, :2], axis=0, weights=weights)
                ax.scatter(
                    *center,
                    marker='x',
                    s=60,
                    color=colors(c - 1),
                    linewidth=2,
                )

        ax.set_title(f'{title} (n = {counts[sample].sum()} reads)')
        ax.set_xlabel(f'PC1 ({pcvar[0]:.1%} variance)')

    y_label = (
        f'PC2 ({pcvar[1]:.1%} variance)'
        if ncomp > 1
        else 'No second variable dimension'
    )
    axes[0].set_ylabel(y_label)
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker='o',
            linestyle='',
            color=colors(c - 1),
            markersize=6,
            label=f'C{c}',
        )
        for c in range(1, k + 1)
    ]
    axes[-1].legend(handles=legend_handles, loc='best', fontsize=8)

    title_prefix = 'Peptide sequence properties'
    fig.suptitle(f'{title_prefix} clusters (Ward, k = {k})')
    fig.text(
        0.5,
        0.015,
        'Point area = fraction of reads (same scale); crosses = read-weighted '
        'cluster centres.\nClustering uses all standardised features; PCA is a '
        'display projection.',
        ha='center',
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 0.94))
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(out / f'pca_comparison.{ext}', dpi=300)
    plt.close(fig)

    # Dendrogram: shows the full hierarchy from which the k clusters were cut.
    fig, ax = plt.subplots(figsize=(12, 5))
    dendrogram(
        Z,
        labels=members.variant_id.to_list(),
        leaf_rotation=90,
        leaf_font_size=7,
        color_threshold=0,
        above_threshold_color='#444444',
        ax=ax,
    )
    ax.set_ylabel('Ward linkage distance')
    ax.set_title('Pooled unique sequences: hierarchical clustering')
    fig.tight_layout()
    fig.savefig(out / 'dendrogram.pdf')
    fig.savefig(out / 'dendrogram.png', dpi=200)
    plt.close(fig)

    # Abundance chart: directly compares each cluster's percentage of reads.
    fig, ax = plt.subplots(figsize=(max(7, k * 1.2), 4))
    x = np.arange(k)
    width = 0.8 / len(sample_names)
    abundance_colors = ['#4575b4', '#66a61e', '#d95f02']
    for sample_index, (sample, title) in enumerate(
        zip(sample_names, sample_labels)
    ):
        offset = (sample_index - (len(sample_names) - 1) / 2) * width
        ax.bar(
            x + offset,
            summary[f'{sample}_fraction'] * 100,
            width,
            label=title,
            color=abundance_colors[sample_index],
        )
    ax.set_xticks(x, [f'C{i}' for i in range(1, k + 1)])
    ax.set_ylabel('Accepted reads (%)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / 'cluster_abundance.pdf')
    fig.savefig(out / 'cluster_abundance.png', dpi=300)
    plt.close(fig)

    # ------------------------------------------------------------------
    # 7. Record reproducibility information and run final consistency checks.
    # ------------------------------------------------------------------
    meta = {
        'features_requested': fs,
        'features_used': active,
        'constant_features_removed': dropped,
        'clusters': k,
        'cluster_selection': a.clusters,
        'silhouette': float(diag.loc[diag.k == k, 'silhouette'].iloc[0]),
        'original_read_totals': original_totals.to_dict(),
        'analysed_read_totals': counts.sum().to_dict(),
        'unique_sequences': len(traits),
        'min_pooled_count': a.min_pooled_count,
        'pca_explained_variance_ratio': pcvar.tolist(),
        'fit_weighting': 'one point per unique sequence',
        'method': 'Ward; Euclidean distance; pooled feature standardisation',
        'sample_order': sample_names,
        'sample_labels': sample_labels,
        'input_files': [str(path) for path in a.inputs],
        # Hashes make it possible to confirm exactly which input files were run.
        'input_sha256': [
            hashlib.sha256(path.read_bytes()).hexdigest()
            for path in a.inputs
        ],
        'versions': {
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'scipy': scipy.__version__,
            'sklearn': sklearn.__version__,
        },
    }
    (out / 'run_metadata.json').write_text(json.dumps(meta, indent=2) + '\n')

    # Every retained read must occur in exactly one cluster.
    for sample in sample_names:
        assert int(summary[f'{sample}_reads'].sum()) == int(counts[sample].sum())

    print(json.dumps(meta, indent=2))
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
