#!/usr/bin/env python3
"""
hotspot_detection.py

Identify variable "hotspot" regions in a position-specific amino-acid
frequency table (the tsv format with columns: position, top_aa, extra_aa,
proportion_breakdown).

Two independent toggles, so you can test what each one is actually buying
you:

  --metric {noise,entropy}
      noise   : 100 - (top_aa %)  (how much of the population disagrees
                with consensus -- ignores *how* the rest is distributed)
      entropy : Shannon entropy over the full residue distribution at that
                position (uses the whole distribution, standard in
                sequence-logo / conservation-scoring tools)

  --detector {mad,plain}
      mad     : robust z-score thresholding -- flag positions where
                (score - median) / (1.4826 * MAD) > k, then merge adjacent
                flagged positions into regions. Robust to the hotspots
                themselves skewing the baseline.
      plain   : naive baseline -- flag positions where score > (mean + k *
                stdev), same merging step. Provided so you can see directly
                whether MAD vs mean/SD changes which positions get flagged
                (the outlier-masking effect described in conversation).

Usage:
    python hotspot_detection.py data1.tsv
    python hotspot_detection.py data1.tsv --metric entropy --k 2.0
    python hotspot_detection.py data1.tsv --detector plain --k 1.5
    python hotspot_detection.py data1.tsv --k-sweep 0.5,1.0,1.5,2.0,2.5,3.0
"""

import argparse
import csv
import math
import re
import statistics
import sys
from pathlib import Path


def parse_table(path):
    """Parse the tsv format into a list of dicts:
    {position, top_aa, residues: {aa: pct, ...}}"""
    rows = []
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            residues = {}
            for chunk in r["proportion_breakdown"].split(","):
                chunk = chunk.strip()
                if not chunk:
                    continue
                aa, pct = chunk.split(":")
                residues[aa.strip()] = float(pct.strip().rstrip("%"))
            rows.append({
                "position": int(r["position"]),
                "top_aa": r["top_aa"].strip(),
                "residues": residues,
            })
    return rows


def score_noise(row):
    return 100.0 - row["residues"][row["top_aa"]]


def score_entropy(row):
    """Shannon entropy in bits over the observed residue distribution."""
    h = 0.0
    for pct in row["residues"].values():
        p = pct / 100.0
        if p > 0:
            h -= p * math.log2(p)
    return h


METRICS = {"noise": score_noise, "entropy": score_entropy}


def mad_zscores(scores):
    med = statistics.median(scores)
    mad = statistics.median(abs(s - med) for s in scores)
    if mad == 0:
        # Degenerate case: fall back to a tiny epsilon so we don't divide by
        # zero when almost every position has an identical score.
        mad = 1e-9
    return [(s - med) / (1.4826 * mad) for s in scores], med, mad


def plain_zscores(scores):
    mean = statistics.mean(scores)
    sd = statistics.stdev(scores) if len(scores) > 1 else 0.0
    if sd == 0:
        sd = 1e-9
    return [(s - mean) / sd for s in scores], mean, sd


DETECTORS = {"mad": mad_zscores, "plain": plain_zscores}


def merge_regions(positions, flagged, max_gap=1):
    """Merge adjacent (within max_gap) flagged positions into contiguous
    regions. positions/flagged are parallel lists in position order."""
    idx = [i for i, f in enumerate(flagged) if f]
    regions = []
    if not idx:
        return regions
    start = prev = idx[0]
    for i in idx[1:]:
        if i - prev > max_gap + 1:
            regions.append((positions[start], positions[prev]))
            start = i
        prev = i
    regions.append((positions[start], positions[prev]))
    return regions


def run(rows, metric, detector, k, max_gap=1):
    positions = [r["position"] for r in rows]
    score_fn = METRICS[metric]
    scores = [score_fn(r) for r in rows]

    z, baseline, spread = DETECTORS[detector](scores)
    flagged = [zi > k for zi in z]
    regions = merge_regions(positions, flagged, max_gap=max_gap)

    return {
        "positions": positions,
        "scores": scores,
        "z": z,
        "baseline": baseline,
        "spread": spread,
        "flagged": flagged,
        "regions": regions,
    }


def print_report(result, rows, metric, detector, k):
    print(f"metric={metric}  detector={detector}  k={k}")
    print(f"baseline ({'median' if detector == 'mad' else 'mean'}) = "
          f"{result['baseline']:.3f}   "
          f"spread ({'MAD' if detector == 'mad' else 'stdev'}) = "
          f"{result['spread']:.3f}\n")

    if not result["regions"]:
        print("No hotspot regions flagged at this threshold.")
        return

    by_pos = {r["position"]: r for r in rows}
    print(f"{len(result['regions'])} hotspot region(s):")
    for start, end in result["regions"]:
        idxs = [i for i, p in enumerate(result["positions"]) if start <= p <= end]
        max_z = max(result["z"][i] for i in idxs)
        vals = [round(result["scores"][i], 1) for i in idxs]
        width = end - start + 1
        shape = "spike" if width == 1 else "block" if width >= 4 else "short run"
        print(f"  positions {start}-{end} ({width} pos, {shape})  "
              f"{metric}={vals}  max_z={max_z:.2f}")


def sweep_k(rows, metric, detector, ks, max_gap=1):
    print(f"k-sweep  metric={metric}  detector={detector}")
    for k in ks:
        result = run(rows, metric, detector, k, max_gap=max_gap)
        n_regions = len(result["regions"])
        n_positions = sum(1 for f in result["flagged"] if f)
        print(f"  k={k:>4}: {n_regions} region(s), {n_positions} flagged position(s)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table", type=Path, help="Path to the tsv frequency table")
    p.add_argument("--metric", choices=METRICS.keys(), default="noise")
    p.add_argument("--detector", choices=DETECTORS.keys(), default="mad")
    p.add_argument("--k", type=float, default=1.5, help="Threshold multiplier (default: 1.5)")
    p.add_argument("--max-gap", type=int, default=1,
                    help="Max gap (in positions) allowed when merging adjacent flagged positions (default: 1)")
    p.add_argument("--k-sweep", type=str, default=None,
                    help="Comma-separated list of k values to sweep instead of a single run, "
                         "e.g. --k-sweep 0.5,1.0,1.5,2.0,2.5,3.0")
    args = p.parse_args()

    if not args.table.exists():
        raise SystemExit(f"File not found: {args.table}")

    rows = parse_table(args.table)

    if args.k_sweep:
        ks = [float(x) for x in args.k_sweep.split(",")]
        sweep_k(rows, args.metric, args.detector, ks, max_gap=args.max_gap)
    else:
        result = run(rows, args.metric, args.detector, args.k, max_gap=args.max_gap)
        print_report(result, rows, args.metric, args.detector, args.k)


if __name__ == "__main__":
    main()