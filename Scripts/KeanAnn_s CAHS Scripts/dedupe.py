#!/usr/bin/env python3
"""
dedupe_sequences.py

Removes duplicate rows from a CSV based on a sequence column, keeping the
first occurrence of each unique sequence. Built for cleaning up scored
mutant libraries (e.g. output from motif1_helix_scorer.py) that may
contain repeats.

USAGE:
    python dedupe_sequences.py input.csv output.csv
    python dedupe_sequences.py input.csv output.csv --column motif1_sequence
    python dedupe_sequences.py input.csv output.csv --case-insensitive

INPUT CSV:
    Any CSV with a sequence-like column. By default the script looks for
    (in order of preference): "sequence", "motif1_sequence", "protein_sequence".
    Use --column to specify a different one explicitly.

OUTPUT CSV:
    Same columns as input, duplicate rows removed. A summary is printed
    (rows in, unique out, duplicates removed).
"""

import sys
import argparse
import pandas as pd

DEFAULT_COLUMN_PRIORITY = ["sequence", "motif1_sequence", "protein_sequence"]


def find_column(columns, requested):
    if requested:
        # Case-insensitive exact match against requested name
        lower_map = {c.lower(): c for c in columns}
        if requested.lower() in lower_map:
            return lower_map[requested.lower()]
        print(f"[error] Column '{requested}' not found. Available columns: {list(columns)}")
        sys.exit(1)

    lower_map = {c.lower(): c for c in columns}
    for candidate in DEFAULT_COLUMN_PRIORITY:
        if candidate in lower_map:
            return lower_map[candidate]

    print(f"[error] Could not auto-detect a sequence column. "
          f"Available columns: {list(columns)}\n"
          f"        Specify one explicitly with --column")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Remove duplicate rows from a CSV based on a sequence column.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("output_file", help="Output CSV file (duplicates removed)")
    parser.add_argument(
        "--column", default=None,
        help="Name of the column to dedupe on (default: auto-detect "
             "'sequence', 'motif1_sequence', or 'protein_sequence')",
    )
    parser.add_argument(
        "--case-insensitive", action="store_true",
        help="Treat sequences differing only in case as duplicates (default: case-sensitive)",
    )
    parser.add_argument(
        "--ignore-whitespace", action="store_true",
        help="Strip whitespace from sequences before comparing (default: off)",
    )
    parser.add_argument(
        "--keep",
        choices=["first", "last"],
        default="first",
        help="Which duplicate copy to keep (default: first)",
    )
    args = parser.parse_args()

    print(f"[read] Reading {args.input_file}...", end="", flush=True)
    try:
        df = pd.read_csv(args.input_file)
    except FileNotFoundError:
        print(f"\n[error] File not found: {args.input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[error] Could not read file: {e}")
        sys.exit(1)
    print(f" {len(df):,} rows")

    if df.empty:
        print("[error] Input CSV is empty")
        sys.exit(1)

    seq_col = find_column(df.columns, args.column)
    print(f"[ok] Deduping on column '{seq_col}'")

    # Build a normalized comparison key without altering the original data
    key = df[seq_col].astype(str)
    if args.ignore_whitespace:
        key = key.str.strip()
    if args.case_insensitive:
        key = key.str.upper()

    n_before = len(df)
    n_blank = int((key.str.len() == 0).sum())

    dedup_mask = ~key.duplicated(keep=args.keep)
    out_df = df[dedup_mask]

    n_after = len(out_df)
    n_removed = n_before - n_after

    print(f"[write] Writing to {args.output_file}...", end="", flush=True)
    try:
        out_df.to_csv(args.output_file, index=False)
    except Exception as e:
        print(f"\n[error] Could not write output: {e}")
        sys.exit(1)
    print(" done")

    print()
    print("=" * 60)
    print("[ok] Deduplication complete!")
    print(f"  Input rows:          {n_before:,}")
    print(f"  Unique rows kept:    {n_after:,}")
    print(f"  Duplicates removed:  {n_removed:,}")
    if n_blank:
        print(f"  Blank/empty sequences: {n_blank:,} (treated as duplicates of each other)")
    print(f"  Output:              {args.output_file}")
    print("=" * 60)


if __name__ == "__main__":
    main()