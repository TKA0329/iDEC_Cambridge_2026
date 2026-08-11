#!/usr/bin/env python3
"""
hotspot_detection.py

Identify variable "hotspot" regions in a set of aligned protein sequences.

Standalone workflow (this is what running the script directly does now):

    aligned FASTA  -->  build_table()  -->  position-specific frequency table
                                             (position, top_aa, extra_aa,
                                             proportion_breakdown)
                                         -->  score + filter (below)

This mirrors the table-building step that filter_sequences_pt2.py performs
inside run_pipeline.py, so hotspot detection no longer depends on a
pipeline run having already produced a pt2_report.tsv -- point it at any
aligned FASTA and it builds the same shape of table itself.

parse_table() (tsv reader) and run() are kept as-is for backward
compatibility: run_pipeline.py imports this module directly and calls
those two functions on its own pt2_report.tsv, so that integration
continues to work unchanged.

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
    python hotspot_detection.py aligned.fasta
    python hotspot_detection.py aligned.fasta --metric entropy --k 2.0
    python hotspot_detection.py aligned.fasta --detector plain --k 1.5
    python hotspot_detection.py aligned.fasta --k-sweep 0.5,1.0,1.5,2.0,2.5,3.0
    python hotspot_detection.py aligned.fasta --save-table freq_table.tsv
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


def parse_fasta(path):
    """Parse a FASTA file into a list of (header, sequence) tuples.

        >1
        seq1
        >2
        seq2

    Sequence lines for a record are concatenated (in case a sequence wraps
    multiple lines). Blank lines are ignored.
    """
    records = []
    header = None
    chunks = []
    with open(path, encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records.append((header, "".join(chunks)))
    return records


def build_table(records, gap_char="-", include_gaps=False):
    """Build the same shape of position-specific frequency table that
    parse_table() reads from tsv, directly from a set of aligned sequences.

    records: list of (header, sequence) tuples, all the same length (i.e.
    already aligned -- this does not perform alignment itself).

    For each column (position), tallies residue composition across all
    sequences and returns a list of dicts:
        {position, top_aa, extra_aa, residues: {aa: pct, ...}}
    (position is 1-indexed, matching the tsv format's convention.)

    By default, gap characters are excluded from the composition at a
    position (percentages are renormalized over the non-gap residues seen
    there) since a gap reflects an alignment artifact, not a "the
    population disagrees with consensus" signal. Pass include_gaps=True to
    instead count the gap character as just another residue.
    """
    if not records:
        raise ValueError("No sequences found -- nothing to build a table from.")

    lengths = {len(seq) for _, seq in records}
    if len(lengths) > 1:
        detail = ", ".join(f"{h}={len(s)}" for h, s in records)
        raise ValueError(
            f"Sequences are not aligned (found {len(lengths)} distinct lengths: "
            f"{sorted(lengths)}). hotspot detection needs an aligned FASTA -- "
            f"align the sequences first. Per-record lengths: {detail}"
        )
    length = lengths.pop()

    rows = []
    for pos in range(length):
        counts = {}
        total = 0
        for _, seq in records:
            aa = seq[pos]
            if aa == gap_char and not include_gaps:
                continue
            counts[aa] = counts.get(aa, 0) + 1
            total += 1
        if total == 0:
            # Every sequence had a gap at this column -- nothing to score.
            continue
        residues = {aa: (count / total) * 100.0 for aa, count in counts.items()}
        ranked = sorted(residues.items(), key=lambda kv: kv[1], reverse=True)
        top_aa = ranked[0][0]
        extra_aa = ranked[1][0] if len(ranked) > 1 else ""
        rows.append({
            "position": pos + 1,
            "top_aa": top_aa,
            "extra_aa": extra_aa,
            "residues": residues,
        })
    return rows


def write_table(rows, path):
    """Write a frequency table (as returned by build_table/parse_table, or
    annotate_table) out as tsv. If rows came from annotate_table() (i.e.
    they carry score/z/hotspot keys), those columns are included too --
    otherwise this writes the same position/top_aa/extra_aa/
    proportion_breakdown format that filter_sequences_pt2.py produces."""
    annotated = rows and "hotspot" in rows[0]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["position", "top_aa", "extra_aa", "proportion_breakdown"]
        if annotated:
            header += ["score", "z", "hotspot"]
        writer.writerow(header)
        for row in rows:
            ranked = sorted(row["residues"].items(), key=lambda kv: kv[1], reverse=True)
            breakdown = ",".join(f"{aa}:{pct:.2f}%" for aa, pct in ranked)
            line = [row["position"], row["top_aa"], row.get("extra_aa", ""), breakdown]
            if annotated:
                line += [f"{row['score']:.4f}", f"{row['z']:.4f}", row["hotspot"]]
            writer.writerow(line)


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

    print(f"{len(result['regions'])} hotspot region(s):")
    for start, end in result["regions"]:
        idxs = [i for i, p in enumerate(result["positions"]) if start <= p <= end]
        max_z = max(result["z"][i] for i in idxs)
        vals = [round(result["scores"][i], 1) for i in idxs]
        width = end - start + 1
        shape = "spike" if width == 1 else "block" if width >= 4 else "short run"

        # Which positions in this span actually crossed k, vs which were
        # only pulled in because they sit within --max-gap of a flagged
        # neighbor. merge_regions() can include the latter in (start, end)
        # without them individually being hotspots.
        hotspot_positions = [result["positions"][i] for i in idxs if result["flagged"][i]]
        gap_filled = [result["positions"][i] for i in idxs if not result["flagged"][i]]

        print(f"  positions {start}-{end} ({width} pos, {shape})  "
              f"{metric}={vals}  max_z={max_z:.2f}")
        print(f"    flagged (z > k): {hotspot_positions}")
        if gap_filled:
            print(f"    pulled in by --max-gap (not individually over k): {gap_filled}")


def annotate_table(rows, result):
    """Merge a frequency table with a run() result into per-position rows
    suitable for write_table(..., annotate=True): adds score/z/hotspot
    columns so every position -- not just merged region bounds -- is
    explicitly marked as a hotspot or not."""
    annotated = []
    for row, score, z, flagged in zip(rows, result["scores"], result["z"], result["flagged"]):
        merged = dict(row)
        merged["score"] = score
        merged["z"] = z
        merged["hotspot"] = flagged
        annotated.append(merged)
    return annotated


def sweep_k(rows, metric, detector, ks, max_gap=1, save_csv=None):
    print(f"k-sweep  metric={metric}  detector={detector}")
    csv_rows = []
    for k in ks:
        result = run(rows, metric, detector, k, max_gap=max_gap)
        n_regions = len(result["regions"])
        n_positions = sum(1 for f in result["flagged"] if f)

        # Build a human-readable region summary, e.g.:
        #   6-10 (pos 6,7,8,9,10), 25 (pos 25), 37-39 (pos 37,38,39)
        region_strs = []
        for start, end in result["regions"]:
            idxs = [i for i, p in enumerate(result["positions"]) if start <= p <= end]
            hotspot_positions = [result["positions"][i] for i in idxs if result["flagged"][i]]
            label = f"{start}" if start == end else f"{start}-{end}"
            region_strs.append(f"{label} (pos {','.join(map(str, hotspot_positions))})")

        print(f"  k={k:>4}: {n_regions} region(s), {n_positions} flagged position(s)")
        if region_strs:
            print(f"        regions: {'; '.join(region_strs)}")

        if save_csv is not None:
            for start, end in result["regions"]:
                idxs = [i for i, p in enumerate(result["positions"]) if start <= p <= end]
                hotspot_positions = [result["positions"][i] for i in idxs if result["flagged"][i]]
                csv_rows.append({
                    "k": k,
                    "region_start": start,
                    "region_end": end,
                    "flagged_positions": ",".join(map(str, hotspot_positions)),
                    "n_flagged": len(hotspot_positions),
                })

    if save_csv is not None:
        with open(save_csv, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["k", "region_start", "region_end", "flagged_positions", "n_flagged"],
                delimiter="\t",
            )
            writer.writeheader()
            for row in csv_rows:
                writer.writerow(row)
        print(f"\nWrote per-k region table ({len(csv_rows)} rows) to {save_csv}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", type=Path, help="Path to an aligned protein FASTA file (all sequences same length)")
    p.add_argument("--metric", choices=METRICS.keys(), default="noise")
    p.add_argument("--detector", choices=DETECTORS.keys(), default="mad")
    p.add_argument("--k", type=float, default=1.5, help="Threshold multiplier (default: 1.5)")
    p.add_argument("--max-gap", type=int, default=1,
                    help="Max gap (in positions) allowed when merging adjacent flagged positions (default: 1)")
    p.add_argument("--k-sweep", type=str, default=None,
                    help="Comma-separated list of k values to sweep instead of a single run, "
                         "e.g. --k-sweep 0.5,1.0,1.5,2.0,2.5,3.0")
    p.add_argument("--gap-char", type=str, default="-",
                    help="Character representing an alignment gap (default: '-')")
    p.add_argument("--include-gaps", action="store_true",
                    help="Count the gap character as a residue when computing per-position "
                         "composition, instead of excluding it and renormalizing over the "
                         "non-gap residues at that column (default: excluded)")
    p.add_argument("--save-table", type=Path, default=None,
                    help="Optional: also write the derived position-frequency table to this "
                         "tsv path (same format filter_sequences_pt2.py's pt2_report.tsv uses)")
    p.add_argument("--save-regions", type=Path, default=None,
                    help="Optional (--k-sweep mode only): write a per-k region/position table "
                         "(k, region_start, region_end, flagged_positions, n_flagged) to this "
                         "tsv path -- this is what actually lists the regions/positions for "
                         "every k in the sweep, since --save-table in sweep mode only writes "
                         "the raw unscored frequency table")
    args = p.parse_args()

    if not args.fasta.exists():
        raise SystemExit(f"File not found: {args.fasta}")

    records = parse_fasta(args.fasta)
    if not records:
        raise SystemExit(f"No sequences found in {args.fasta}")

    try:
        rows = build_table(records, gap_char=args.gap_char, include_gaps=args.include_gaps)
    except ValueError as exc:
        raise SystemExit(str(exc))

    if args.save_table and args.k_sweep:
        # No single result to annotate against when sweeping k -- write the
        # raw (unscored) frequency table instead.
        write_table(rows, args.save_table)
        print(f"Wrote position-frequency table ({len(rows)} positions) to {args.save_table}\n")

    if args.k_sweep:
        ks = [float(x) for x in args.k_sweep.split(",")]
        sweep_k(rows, args.metric, args.detector, ks, max_gap=args.max_gap,
                save_csv=args.save_regions)
    else:
        result = run(rows, args.metric, args.detector, args.k, max_gap=args.max_gap)
        print_report(result, rows, args.metric, args.detector, args.k)
        if args.save_table:
            write_table(annotate_table(rows, result), args.save_table)
            n_hotspot = sum(result["flagged"])
            print(f"\nWrote per-position table with score/z/hotspot flag "
                  f"({len(rows)} positions, {n_hotspot} flagged) to {args.save_table}")


if __name__ == "__main__":
    main()