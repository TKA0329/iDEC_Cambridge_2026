#!/usr/bin/env python3
"""
position_frequency_analyzer.py

Computes the numeric version of your WebLogo: for each position in
motif1_sequence across your whole shortlisted CSV, counts how often each
amino acid appears, and reports the smallest set of residues needed to
cover a given fraction (default 90%) of survivors at that position.

This is meant to run on your ROUND 2 shortlist (the CSV with a
motif1_sequence column, e.g. output of motif1_filter_pipeline.py), since
that's the data that should drive your actual degenerate-codon choices per
hotspot position.

USAGE:
    python position_frequency_analyzer.py shortlist.csv position_frequencies.csv

OUTPUT CSV columns:
    position, amino_acid, count, percentage, cumulative_percentage

Also prints, per position, a short summary: how many distinct residues
appear, and the minimal residue set needed to cover >= COVERAGE_THRESHOLD
of all survivors at that position (useful for picking a degenerate codon
narrower than full NNS where the data supports it).
"""

import sys
import csv
from collections import Counter

COVERAGE_THRESHOLD = 0.90  # report the smallest residue set covering >= this fraction


def main():
    if len(sys.argv) != 3:
        print("Usage: python position_frequency_analyzer.py shortlist.csv position_frequencies.csv")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        if "motif1_sequence" not in reader.fieldnames:
            print("ERROR: input CSV must have a 'motif1_sequence' column "
                  "(this should be a scored shortlist, not the raw mutator output).")
            sys.exit(1)

        motif_length = None
        position_counters = None
        n_sequences = 0

        for row in reader:
            seq = row.get("motif1_sequence", "")
            if not seq:
                continue
            if motif_length is None:
                motif_length = len(seq)
                position_counters = [Counter() for _ in range(motif_length)]
            elif len(seq) != motif_length:
                # skip malformed/length-mismatched rows rather than crash
                continue

            for i, aa in enumerate(seq):
                position_counters[i][aa] += 1
            n_sequences += 1

    if n_sequences == 0:
        print("No valid motif1_sequence rows found.")
        sys.exit(1)

    print(f"Analyzed {n_sequences:,} sequences, motif length {motif_length}\n")

    out_rows = []
    for pos_idx, counter in enumerate(position_counters):
        position = pos_idx + 1  # 1-indexed to match your WebLogo
        total = sum(counter.values())
        ranked = counter.most_common()

        cumulative = 0
        cumulative_pct = 0.0
        min_set = []
        for aa, count in ranked:
            cumulative += count
            cumulative_pct = cumulative / total
            min_set.append(aa)
            if cumulative_pct >= COVERAGE_THRESHOLD:
                break

        print(f"Position {position:2d}: {len(ranked)} distinct residues observed. "
              f"Top set covering >={int(COVERAGE_THRESHOLD*100)}%: {','.join(min_set)} "
              f"({cumulative_pct*100:.1f}% actual)")

        running_cum = 0
        for aa, count in ranked:
            running_cum += count
            out_rows.append({
                "position": position,
                "amino_acid": aa,
                "count": count,
                "percentage": round(100 * count / total, 2),
                "cumulative_percentage": round(100 * running_cum / total, 2),
            })

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["position", "amino_acid", "count", "percentage", "cumulative_percentage"])
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"\nFull per-position frequency table -> {out_path}")


if __name__ == "__main__":
    main()