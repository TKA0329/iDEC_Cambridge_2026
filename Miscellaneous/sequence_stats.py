#!/usr/bin/env python3
"""
Length-based statistics for protein sequence CSV files.

Usage:
  python sequence_stats.py input.csv [output_dir]
"""

import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = [
    "sequence",
]

OPTIONAL_COLUMNS = [
    "gravy",
    "aromaticity",
    "instability_index",
    "charge_at_pH7",
    "positive_res_RK",
    "negative_res_DE",
    "helix_fraction",
    "sheet_fraction",
    "turn_fraction",
    "fcr",
    "ncpr",
    "cider_kappa",
]

SUMMARY_COLUMNS = [
    "gravy",
    "aromaticity",
    "instability_index",
    "charge_at_pH7",
    "charge_sum",
    "helix_sheet_sum",
    "turn_fraction",
    "fcr",
    "ncpr",
    "cider_kappa",
]

DERIVED_COLUMN_SOURCES = {
    "charge_sum": ("positive_res_RK", "negative_res_DE"),
    "helix_sheet_sum": ("helix_fraction", "sheet_fraction"),
}


def add_derived_columns(df):
    df = df.copy()
    skipped = {}

    if all(col in df.columns for col in DERIVED_COLUMN_SOURCES["charge_sum"]):
        df["charge_sum"] = df["positive_res_RK"].astype(float) + df["negative_res_DE"].astype(float)
    else:
        skipped["charge_sum"] = [
            col for col in DERIVED_COLUMN_SOURCES["charge_sum"] if col not in df.columns
        ]

    if all(col in df.columns for col in DERIVED_COLUMN_SOURCES["helix_sheet_sum"]):
        df["helix_sheet_sum"] = df["helix_fraction"].astype(float) + df["sheet_fraction"].astype(float)
    else:
        skipped["helix_sheet_sum"] = [
            col for col in DERIVED_COLUMN_SOURCES["helix_sheet_sum"] if col not in df.columns
        ]

    return df, skipped


def validate_columns(df):
    missing_required = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_required:
        raise ValueError(f"Missing required columns: {', '.join(missing_required)}")
    return [c for c in OPTIONAL_COLUMNS if c not in df.columns]


def get_available_summary_columns(df):
    return [col for col in SUMMARY_COLUMNS if col in df.columns]


def build_length_summary(df, length_col):
    row_counts = df.groupby(length_col).size().rename("row_count").reset_index()
    stats = get_available_summary_columns(df)
    if not stats:
        return row_counts

    agg_map = {col: ["min", "max", "mean"] for col in stats}
    grouped = df.groupby(length_col).agg(agg_map)
    grouped.columns = [f"{col}_{stat}" for col, stat in grouped.columns]
    grouped = grouped.reset_index()
    return row_counts.merge(grouped, on=length_col, how="left")


def print_summary_table(summary_df):
    print("\nLength summary stats:")
    print(summary_df.to_string(index=False, float_format="{:.4f}".format))
    print()


def write_summary_csv(summary_df, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, index=False)


def main():
    if len(sys.argv) < 2:
        print("Usage: python sequence_stats.py <input.csv> [output_dir]")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("stats_output")

    if not input_file.exists():
        print(f"Input file not found: {input_file}")
        sys.exit(1)

    print(f"Reading {input_file}")
    df = pd.read_csv(input_file, comment="#")
    if df.empty:
        print("Input CSV is empty")
        sys.exit(1)

    try:
        missing_optional = validate_columns(df)
    except ValueError as exc:
        print(exc)
        sys.exit(1)

    if missing_optional:
        print(f"Warning: missing optional columns will be skipped: {', '.join(missing_optional)}")

    df, skipped_derived = add_derived_columns(df)
    for derived_col, missing_sources in skipped_derived.items():
        print(
            f"Warning: could not compute '{derived_col}' because these source columns are missing: "
            f"{', '.join(missing_sources)}"
        )

    length_col = "length" if "length" in df.columns else df.columns[0]
    if length_col != "length":
        print(f"Warning: 'length' column was not found. Using first column '{length_col}' as length.")

    available_summary_columns = get_available_summary_columns(df)
    if available_summary_columns:
        print(f"Summary metrics included: {', '.join(available_summary_columns)}")
    else:
        print("Warning: no metric columns were available for aggregation. Summary will include row counts only.")

    summary_df = build_length_summary(df, length_col)
    summary_csv_path = output_dir / "length_summary_stats.csv"
    write_summary_csv(summary_df, summary_csv_path)

    print_summary_table(summary_df)
    print(f"Length summary saved to: {summary_csv_path}")


if __name__ == "__main__":
    main()
