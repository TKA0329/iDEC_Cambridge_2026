#!/usr/bin/env python3
"""
length_counter.py

Reads a protein-analyzer CSV (with '#'-prefixed metadata header lines) and:
  1. Counts how many sequences exist at each length.
  2. Writes a length,count CSV.
  3. Optionally reports how many sequences fall at/under given length thresholds.

Usage:
  python length_counter.py input.csv output.csv
  python length_counter.py input.csv output.csv --thresholds 24,50,100,150
  python length_counter.py input.csv output.csv --thresholds 24,50,100,150 --thresholds-out under_length_summary.csv

The input file's leading lines starting with '#' are treated as metadata
and skipped automatically; the first non-'#' line is used as the header.
"""

import argparse
import csv
import sys
from collections import Counter


def read_rows(input_path):
    """Yield DictReader rows from the CSV, skipping leading '#' metadata lines."""
    with open(input_path, newline="", encoding="utf-8") as f:
        # Skip metadata/comment lines until we hit the real header.
        # readline() (unlike the "for line in f" iterator) keeps file
        # position tracking valid, so no seek() gymnastics are needed.
        while True:
            line = f.readline()
            if line == "":
                raise ValueError("No header row found (file was all comments or empty).")
            if not line.startswith("#"):
                header_line = line
                break

        # Feed the header line plus the rest of the file into DictReader
        def line_stream():
            yield header_line
            for remaining_line in f:
                yield remaining_line

        reader = csv.DictReader(line_stream())
        for row in reader:
            yield row


def count_lengths(input_path, length_field="length"):
    counts = Counter()
    total_rows = 0
    skipped = 0

    for row in read_rows(input_path):
        total_rows += 1
        raw = row.get(length_field, "")
        if raw is None or raw == "":
            skipped += 1
            continue
        try:
            length = int(float(raw))  # handles "24" or "24.0"
        except ValueError:
            skipped += 1
            continue
        counts[length] += 1

    return counts, total_rows, skipped


def write_length_counts(counts, output_path):
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["length", "count"])
        for length in sorted(counts):
            writer.writerow([length, counts[length]])


def write_threshold_summary(counts, thresholds, output_path):
    total = sum(counts.values())
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["max_length", "count_at_or_under", "fraction_of_total"])
        for t in sorted(thresholds):
            count_under = sum(c for length, c in counts.items() if length <= t)
            frac = count_under / total if total else 0
            writer.writerow([t, count_under, f"{frac:.4f}"])


def parse_thresholds(arg):
    if not arg:
        return []
    return [int(x.strip()) for x in arg.split(",") if x.strip()]


def main():
    parser = argparse.ArgumentParser(description="Count protein sequences by length.")
    parser.add_argument("input_csv", help="Path to the input CSV (protein analyzer output).")
    parser.add_argument("output_csv", help="Path to write length,count CSV.")
    parser.add_argument(
        "--thresholds",
        default="",
        help="Comma-separated list of lengths to report cumulative counts under, e.g. '24,50,100,150'.",
    )
    parser.add_argument(
        "--thresholds-out",
        default=None,
        help="Path to write the threshold summary CSV (default: <output_csv> with '_thresholds' suffix).",
    )
    parser.add_argument(
        "--length-field",
        default="length",
        help="Name of the length column in the input CSV (default: 'length').",
    )

    args = parser.parse_args()

    print(f"Reading {args.input_csv} ...")
    counts, total_rows, skipped = count_lengths(args.input_csv, args.length_field)
    print(f"Processed {total_rows:,} rows ({skipped:,} skipped due to missing/invalid length).")

    write_length_counts(counts, args.output_csv)
    print(f"Wrote length/count table to {args.output_csv} ({len(counts)} distinct lengths).")

    thresholds = parse_thresholds(args.thresholds)
    if thresholds:
        out_path = args.thresholds_out
        if out_path is None:
            if args.output_csv.lower().endswith(".csv"):
                out_path = args.output_csv[:-4] + "_thresholds.csv"
            else:
                out_path = args.output_csv + "_thresholds.csv"
        write_threshold_summary(counts, thresholds, out_path)
        print(f"Wrote threshold summary to {out_path}.")


if __name__ == "__main__":
    main()