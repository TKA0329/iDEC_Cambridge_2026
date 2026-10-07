#!/usr/bin/env python3
"""Plot cluster abundance (A) and weighted trait means (B) for 2-3 samples.

Reads cluster_counts.csv, cluster_traits.csv, and members.csv from ONE clustering
run. Does not recluster, fit PCA, calculate p-values, or infer replication.
Install: python -m pip install numpy pandas matplotlib
Run: python plot_cluster_changes.py results --outdir panels --title SSB
"""
import argparse
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = ['#9fc8c8', '#65a6a6', '#176f6f']
SAMPLE_ORDER = ['unselected', 'intermediate', 'selected']
LABELS = {'gravy':'GRAVY', 'charge_at_pH7':'Estimated net charge at pH 7',
          'instability_index':'Instability index', 'turn_fraction':'Turn fraction',
          'helix_sheet_sum':'Helix + sheet fraction', 'pred_cv_JperK':'Predicted heat capacity (J/K)', 'pred_rog_A':'Predicted radius of gyration (Å)', 'pred_tau_fs':'Predicted end-to-end decorrelation time (fs)'}


def read_outputs(directory):
    """Validate the summaries and return their two or three sample names."""
    counts = pd.read_csv(directory/'cluster_counts.csv', dtype={'cluster':str})
    traits = pd.read_csv(directory/'cluster_traits.csv', dtype={'cluster':str})
    members = pd.read_csv(directory/'members.csv', dtype={'cluster':str})
    samples = [s for s in SAMPLE_ORDER if s+'_reads' in counts]
    if len(samples) not in (2, 3) or samples[0] != 'unselected' or samples[-1] != 'selected':
        raise ValueError('Counts must contain unselected/selected and may contain intermediate')
    needed = {'cluster','feature',*[s+'_read_weighted_mean' for s in samples]}
    if not needed.issubset(traits): raise ValueError(f'Missing trait columns: {needed-set(traits)}')
    if counts.empty or counts.cluster.isna().any() or counts.cluster.duplicated().any():
        raise ValueError('Counts must contain one nonempty row per cluster')
    if traits[['cluster','feature']].isna().any().any() or traits.duplicated(['cluster','feature']).any():
        raise ValueError('Traits must contain one row per cluster/property pair')
    if set(counts.cluster)!=set(traits.cluster): raise ValueError('Cluster IDs disagree between files')
    if members.empty or members.cluster.isna().any(): raise ValueError('Members must have cluster IDs')
    if set(counts.cluster)!=set(members.cluster): raise ValueError('Member cluster IDs disagree with summaries')
    counts = counts.set_index('cluster')
    for sample in samples:
        col = sample+'_reads'; counts[col] = pd.to_numeric(counts[col], errors='raise')
        members[col] = pd.to_numeric(members[col], errors='raise')
        x = counts[col].to_numpy()
        member_x = members[col].to_numpy()
        if not np.isfinite(x).all() or (x<0).any() or (x%1!=0).any() or x.sum()==0:
            raise ValueError(f'{col}: expected nonnegative integer counts and positive total')
        if not np.isfinite(member_x).all() or (member_x<0).any() or (member_x%1!=0).any():
            raise ValueError(f'members/{col}: expected nonnegative integer counts')
        member_totals = members.groupby('cluster')[col].sum().reindex(counts.index)
        if member_totals.isna().any() or not np.array_equal(member_totals.to_numpy(), x):
            raise ValueError(f'{col}: member counts disagree with cluster counts')
        # Recompute proportions from ALL cluster counts. No renormalisation to
        # a selected subset and no use of unique-sequence counts as abundance.
        fraction = counts[col]/counts[col].sum()
        if sample+'_fraction' in counts and not np.allclose(counts[sample+'_fraction'], fraction):
            raise ValueError('Stored proportions disagree with read counts')
        counts[sample+'_fraction'] = fraction
    return counts, traits, members, samples


def save(fig, directory, name):
    """Raster preview plus vector formats suitable for report assembly."""
    for ext in ['png','pdf','svg']:
        fig.savefig(directory/f'{name}.{ext}', dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results_dir', type=Path)
    parser.add_argument('--outdir',type=Path,default=Path('cluster_panels'))
    parser.add_argument('--title',default='Tag library')
    parser.add_argument('--before-label',default='Untreated')
    parser.add_argument('--middle-label',default='Intermediate')
    parser.add_argument('--after-label',default='Selected')
    parser.add_argument('--features',help='Optional comma-separated subset, default all features in traits file')
    # Retained so existing commands using --show-sd continue to work. Weighted
    # SD error bars are now always shown.
    parser.add_argument('--show-sd',action='store_true',help=argparse.SUPPRESS)
    args = parser.parse_args()
    counts,traits,members,samples=read_outputs(args.results_dir)
    clusters=counts.index.tolist()
    features=list(dict.fromkeys(traits.feature)) if not args.features else [f.strip() for f in args.features.split(',')]
    if not features or not all(features) or len(set(features))!=len(features): raise ValueError('Invalid feature list')
    if not set(features).issubset(traits.feature): raise ValueError('Requested feature absent from traits')
    if not set(features).issubset(members.columns): raise ValueError('Requested feature absent from members')
    for feature in features:
        members[feature]=pd.to_numeric(members[feature],errors='raise')
        if not np.isfinite(members[feature]).all(): raise ValueError(f'{feature}: nonfinite member values')
    args.outdir.mkdir(parents=True,exist_ok=True)
    names=(
        [args.before_label,args.after_label]
        if len(samples)==2
        else [args.before_label,args.middle_label,args.after_label]
    )
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})

    # PANEL A: fraction of the whole sequenced sample belonging to each cluster.
    # Bars use different sample-specific denominators.
    fig,ax=plt.subplots(figsize=(max(7,len(clusters)*1.0),4.8),layout='constrained')
    x=np.arange(len(clusters)); width=.8/len(samples)
    for j,sample in enumerate(samples):
        vals=100*counts[sample+'_fraction'].to_numpy()
        offset=(j-(len(samples)-1)/2)*width
        bars=ax.bar(x+offset,vals,width,color=COLORS[j],
                    label=f'{names[j]} (n = {int(counts[sample+"_reads"].sum())} reads)')
        ax.bar_label(bars,labels=[f'{v:.1f}%' for v in vals],padding=3,fontsize=9)
    ax.set_xticks(x,['C'+c for c in clusters]);ax.set_ylim(0,115)
    ax.set_yticks([0,20,40,60,80,100]);ax.set_ylabel('Accepted reads (%)')
    ax.set_title(f'A. {args.title}: cluster representation');ax.legend(loc='upper left',fontsize=9)
    ax.spines[['top','right']].set_visible(False)
    save(fig,args.outdir,'panel_A_cluster_abundance')
    a=counts.copy();a['percentage_point_change']=100*(a.selected_fraction-a.unselected_fraction)
    for previous,current in zip(samples,samples[1:]):
        a[current+'_minus_'+previous+'_percentage_points']=100*(
            a[current+'_fraction']-a[previous+'_fraction']
        )
    a.to_csv(args.outdir/'panel_A_values.csv')

    # PANEL B: conditional mean and population SD within each cluster, both
    # calculated directly from members.csv and weighted by sample read counts.
    # No reads => undefined statistics, displayed as absent rather than zero.
    ncols=min(2,len(features));nrows=math.ceil(len(features)/ncols)
    fig,axes=plt.subplots(nrows,ncols,figsize=(11,3.2*nrows+1),squeeze=False)
    exports=[]
    for ax,feature in zip(axes.flat,features):
        sub=traits[traits.feature==feature].set_index('cluster')
        if set(sub.index)!=set(clusters): raise ValueError(f'{feature}: missing clusters')
        sub=sub.reindex(clusters)
        offsets=np.linspace(-.18,.18,len(samples))
        for i,c in enumerate(clusters):
            cluster_members=members[members.cluster==c]
            values=cluster_members[feature].to_numpy()
            means=[]
            row={'cluster':c,'feature':feature}
            for j,sample in enumerate(samples):
                position=i+offsets[j]
                weights=cluster_members[sample+'_reads'].to_numpy()
                present=weights.sum()>0
                if present:
                    mu=float(np.average(values,weights=weights))
                    sd=float(np.sqrt(np.average((values-mu)**2,weights=weights)))
                    stored_mu=float(sub.loc[c,sample+'_read_weighted_mean'])
                    if not np.isclose(mu,stored_mu,rtol=1e-10,atol=1e-12):
                        raise ValueError(f'{c}/{feature}/{sample}: member mean disagrees with traits')
                else:
                    mu=sd=np.nan
                means.append(mu)
                row[sample+'_read_weighted_mean']=mu
                row[sample+'_read_weighted_sd']=sd
                if present:
                    ax.scatter(position,mu,s=45,color=COLORS[j],zorder=3)
                else:
                    ax.text(position,.03,'absent',rotation=90,
                            transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color=COLORS[j])
            finite=np.isfinite(means)
            if finite.sum()>1:
                ax.plot(i+offsets[finite],np.asarray(means)[finite],color='#999999',linewidth=1,zorder=1)
            for previous_index,(previous,current) in enumerate(zip(samples,samples[1:])):
                row[current+'_minus_'+previous]=means[previous_index+1]-means[previous_index]
            row['selected_minus_unselected']=means[-1]-means[0];exports.append(row)
        ax.set_xticks(x,['C'+c for c in clusters]);ax.set_xlim(-.5,len(clusters)-.5)
        ax.set_title(LABELS.get(feature,feature));ax.set_ylabel('Read-weighted mean')
        ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
    for ax in axes.flat[len(features):]:ax.remove()
    fig.suptitle(f'B. {args.title}: traits within each fixed cluster',y=.99,fontsize=14)
    fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color=COLORS[j],label=names[j]) for j in range(len(samples))],
               loc='upper center',bbox_to_anchor=(.5,.96),ncol=len(samples),frameon=False)
    note='Each point is a within-cluster read-weighted mean; connecting lines show changes across samples.'
    fig.text(.5,.015,note,ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.07,1,.91))
    save(fig,args.outdir,'panel_B_weighted_traits')
    pd.DataFrame(exports).to_csv(args.outdir/'panel_B_values.csv',index=False)
    print(f'Saved panels A/B as PNG, PDF and SVG plus plotted-value CSVs to {args.outdir}')

if __name__=='__main__':main()
