#!/usr/bin/env python3
"""
extract_glycine_motifs.py

Pulls out just the rows from a scored CSV (output of motif1_filter_pipeline.py
or motif1_helix_scorer.py) where motif1_sequence contains at least one
Glycine, and adds columns describing WHERE in the motif that Gly sits --
so you can actually check the "most survivors have Gly in the middle"
pattern systematically instead of eyeballing it.

Adds:
    glycine_positions   -> 1-indexed position(s) of G within motif1_sequence,
                           comma-separated if more than one (e.g. "5" or "5,12")
    glycine_region       -> front / middle / back, based on which third of
                           motif1's length each G falls into. If a sequence
                           has multiple Gly in different regions, all
                           regions are listed (e.g. "front,middle")

USAGE:
    python extract_glycine_motifs.py scored_input.csv glycine_only_output.csv
"""

import sys
import csv


def classify_region(position: int, motif_length: int) -> str:
    """1-indexed position -> 'front', 'middle', or 'back' third of the motif."""
    third = motif_length / 3.0
    if position <= third:
        return "front"
    elif position <= 2 * third:
        return "middle"
    else:
        return "back"


def find_glycine_info(motif1_seq: str):
    positions = [i + 1 for i, aa in enumerate(motif1_seq) if aa == "G"]  # 1-indexed
    regions = []
    for pos in positions:
        region = classify_region(pos, len(motif1_seq))
        if region not in regions:
            regions.append(region)
    return positions, regions


def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_glycine_motifs.py scored_input.csv glycine_only_output.csv")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        if "motif1_sequence" not in reader.fieldnames:
            print("ERROR: input CSV must have a 'motif1_sequence' column "
                  "(this should be the output of motif1_filter_pipeline.py "
                  "or motif1_helix_scorer.py, not the raw mutator output).")
            sys.exit(1)

        rows_out = []
        region_counts = {"front": 0, "middle": 0, "back": 0}

        for row in reader:
            motif1 = row.get("motif1_sequence", "")
            if "G" not in motif1:
                continue

            positions, regions = find_glycine_info(motif1)
            for region in regions:
                region_counts[region] += 1

            new_row = dict(row)
            new_row["glycine_positions"] = ",".join(str(p) for p in positions)
            new_row["glycine_region"] = ",".join(regions)
            rows_out.append(new_row)

    if not rows_out:
        print("No glycine-containing motif1 sequences found in input.")
        sys.exit(0)

    fieldnames = list(rows_out[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Found {len(rows_out):,} glycine-containing motif1 sequences -> {out_path}")
    print("Region breakdown (a sequence with Gly in multiple regions counts in each):")
    for region, count in region_counts.items():
        print(f"  {region}: {count:,}")


if __name__ == "__main__":
    main()