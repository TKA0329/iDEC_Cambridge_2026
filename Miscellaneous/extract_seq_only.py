#!/usr/bin/env python3
"""
extract_sequences.py

Pulls just the 'sequence' column out of a batch_mutator.py output CSV
and writes a new CSV containing only that single 'sequence' column.

USAGE:
    python extract_sequences.py input.csv output.csv
"""

import sys
import csv


def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_sequences.py input.csv output.csv")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        if "sequence" not in reader.fieldnames:
            print(f"ERROR: no 'sequence' column found. Columns present: {reader.fieldnames}")
            sys.exit(1)

        sequences = []
        n_skipped_empty = 0

        for row in reader:
            seq = (row.get("sequence") or "").strip()
            if not seq:
                n_skipped_empty += 1
                continue
            sequences.append(seq)

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"])
        for seq in sequences:
            writer.writerow([seq])

    print(f"Extracted {len(sequences):,} sequences -> {out_path}")
    if n_skipped_empty:
        print(f"  (skipped {n_skipped_empty:,} rows with blank sequence)")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()