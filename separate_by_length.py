#!/usr/bin/env python3
"""
Separate a CSV file by sequence length.
Creates one CSV per length with headers included.

Run: python separate_by_length.py input.csv [output_dir]
"""

import sys
import os
import pandas as pd
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print("Usage: python separate_by_length.py <input.csv> [output_dir]")
        print()
        print("Example:")
        print("  python separate_by_length.py results.csv length_splits/")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "length_splits"

    print(f"Reading {input_file}...")
    try:
        # comment='#' tells pandas to skip any leading '#'-prefixed
        # metadata lines (like the batch_analyzer.py header block)
        # instead of trying to parse them as data/header rows.
        df = pd.read_csv(input_file, comment="#")
    except FileNotFoundError:
        print(f"File not found: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    if df.empty:
        print("Input CSV is empty")
        sys.exit(1)

    # Get sequence column (first column)
    seq_col = df.columns[0]
    print(f"Found {len(df):,} sequences in column '{seq_col}'")

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}/")

    # Calculate sequence lengths
    df["_seq_length"] = df[seq_col].astype(str).str.len()

    # Group by length
    lengths = sorted(df["_seq_length"].unique())
    print(f"Found {len(lengths)} unique lengths: {min(lengths)}-{max(lengths)} amino acids")
    print()

    # Save each length to separate file
    for length in lengths:
        subset = df[df["_seq_length"] == length].drop(columns=["_seq_length"])
        output_file = os.path.join(output_dir, f"length_{length:02d}.csv")
        subset.to_csv(output_file, index=False)
        print(f"  length_{length:02d}.csv: {len(subset):,} sequences")

    print()
    print("=" * 70)
    print(f"Separation complete!")
    print(f"  Total sequences: {len(df):,}")
    print(f"  Unique lengths: {len(lengths)}")
    print(f"  Output directory: {output_dir}/")
    print("=" * 70)


if __name__ == "__main__":
    main()