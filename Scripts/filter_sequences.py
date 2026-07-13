#!/usr/bin/env python3
"""
Interactive CSV filter for protein sequence results.

This script adds derived summary columns when possible, prints per-length stats
to the terminal, writes a length-by-length summary CSV, and generates a filtered
output CSV with filter parameters embedded as a metadata section.

Run:
  python filter_sequences.py input.csv

Optional:
  python filter_sequences.py input.csv output_dir
"""

import re
import sys
from pathlib import Path

import pandas as pd
from sequence_stats import add_derived_columns, build_length_summary, print_summary_table, validate_columns


def prompt_float(prompt_text, default=None, allow_empty=False):
    while True:
        prompt = f"{prompt_text}"
        if default is not None:
            prompt += f" [{default}]"
        prompt += ": "
        value = input(prompt).strip()
        if value == "" and default is not None:
            return float(default)
        if value == "" and allow_empty:
            return None
        try:
            return float(value)
        except ValueError:
            print("Please enter a valid number.")


def prompt_range(prompt_text, default_min=None, default_max=None, allow_empty=False):
    pattern = re.compile(r"^\s*([+-]?\d+(?:\.\d+)?)\s*-\s*([+-]?\d+(?:\.\d+)?)\s*$")
    while True:
        default_text = None
        if default_min is not None and default_max is not None:
            default_text = f"{default_min}-{default_max}"
        prompt = f"{prompt_text}"
        if default_text:
            prompt += f" [{default_text}]"
        prompt += ": "
        value = input(prompt).strip()
        if value == "" and default_text is not None:
            value = default_text
        if value == "" and allow_empty:
            return None, None
        match = pattern.match(value)
        if not match:
            print("Enter a numeric range using a hyphen, e.g. -1.0--0.5 or 6.0-7.0.")
            continue
        try:
            low = float(match.group(1))
            high = float(match.group(2))
        except ValueError:
            print("Please enter valid numeric values.")
            continue
        if low > high:
            print("Minimum cannot be greater than maximum.")
            continue
        return low, high


def write_output_with_metadata(output_path, filtered_df, metadata):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        f.write("# Filtered protein sequence results\n")
        for key, value in metadata.items():
            f.write(f"# {key}: {value}\n")
        f.write("#\n")
        filtered_df.to_csv(f, index=False)


def apply_upper_bound_filter(df, remaining_mask, column, value, report_name):
    next_mask = remaining_mask & (df[column].astype(float) <= value)
    removed = int(remaining_mask.sum() - next_mask.sum())
    return next_mask, {
        report_name: {
            "value": value,
            "removed": removed,
            "remaining_after": int(next_mask.sum()),
        }
    }


def apply_lower_bound_filter(df, remaining_mask, column, value, report_name):
    next_mask = remaining_mask & (df[column].astype(float) > value)
    removed = int(remaining_mask.sum() - next_mask.sum())
    return next_mask, {
        report_name: {
            "value": value,
            "removed": removed,
            "remaining_after": int(next_mask.sum()),
        }
    }


def apply_range_filter(df, remaining_mask, column, low, high, report_name):
    next_mask = remaining_mask & (df[column].astype(float) >= low)
    next_mask &= df[column].astype(float) <= high
    removed = int(remaining_mask.sum() - next_mask.sum())
    return next_mask, {
        report_name: {
            "value": f"{low}-{high}",
            "removed": removed,
            "remaining_after": int(next_mask.sum()),
        }
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python filter_sequences.py <input.csv> [output_dir]")
        print()
        print("Example:")
        print("  python filter_sequences.py results.csv filtered_output")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("filtered_output")

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

    summary_df = build_length_summary(df, length_col)
    summary_csv_path = output_dir / "length_summary_stats.csv"
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(summary_csv_path, index=False)

    print_summary_table(summary_df)
    print(f"Length summary saved to: {summary_csv_path}\n")

    available_filter_descriptions = []
    skipped_filters = []

    if "gravy" in df.columns:
        available_filter_descriptions.append("gravy: keep rows with gravy <= value")
    else:
        skipped_filters.append("gravy")
    if "instability_index" in df.columns:
        available_filter_descriptions.append("instability_index: keep rows with instability_index <= value")
    else:
        skipped_filters.append("instability_index")
    if "helix_sheet_sum" in df.columns:
        available_filter_descriptions.append(
            "helix_sheet_sum: keep rows with helix_fraction + sheet_fraction <= value"
        )
    else:
        skipped_filters.append("helix_sheet_sum")
    if "turn_fraction" in df.columns:
        available_filter_descriptions.append("turn_fraction: keep rows with turn_fraction > value")
    else:
        skipped_filters.append("turn_fraction")
    if "charge_sum" in df.columns:
        available_filter_descriptions.append(
            "charge_sum: keep rows with positive_res_RK + negative_res_DE > value"
        )
    else:
        skipped_filters.append("charge_sum")
    if "charge_at_pH7" in df.columns:
        available_filter_descriptions.append("charge_at_pH7: keep rows within the specified range")
    else:
        skipped_filters.append("charge_at_pH7")
    if "fcr" in df.columns:
        available_filter_descriptions.append("fcr: keep rows with FCR <= value")
    else:
        skipped_filters.append("fcr")
    if "ncpr" in df.columns:
        available_filter_descriptions.append("ncpr: keep rows with NCPR within range")
    else:
        skipped_filters.append("ncpr")
    if "cider_kappa" in df.columns:
        available_filter_descriptions.append("cider_kappa: keep rows with CIDER kappa <= value")
    else:
        skipped_filters.append("cider_kappa")

    print("Enter filter values. Press Enter to skip any available filter.")
    for description in available_filter_descriptions:
        print(f"  - {description}")
    if skipped_filters:
        print(f"  - unavailable because columns are missing: {', '.join(skipped_filters)}")
    print("  - leave blank to skip a filter")
    print()

    gravy_max = prompt_float("gravy <=", allow_empty=True) if "gravy" in df.columns else None
    instability_max = (
        prompt_float("instability_index <=", allow_empty=True)
        if "instability_index" in df.columns
        else None
    )
    helix_sheet_max = (
        prompt_float("helix_fraction + sheet_fraction <=", allow_empty=True)
        if "helix_sheet_sum" in df.columns
        else None
    )
    turn_fraction_min = (
        prompt_float("turn_fraction >", allow_empty=True) if "turn_fraction" in df.columns else None
    )
    charge_sum_min = (
        prompt_float("positive_res_RK + negative_res_DE >", allow_empty=True)
        if "charge_sum" in df.columns
        else None
    )
    if "charge_at_pH7" in df.columns:
        charge_pH7_min, charge_pH7_max = prompt_range(
            "charge_at_pH7 range (min-max)", allow_empty=True
        )
    else:
        charge_pH7_min, charge_pH7_max = None, None
    fcr_max = prompt_float("fcr <=", allow_empty=True) if "fcr" in df.columns else None
    if "ncpr" in df.columns:
        ncpr_min, ncpr_max = prompt_range("ncpr range (min-max)", allow_empty=True)
    else:
        ncpr_min, ncpr_max = None, None
    cider_kappa_max = (
        prompt_float("cider_kappa <=", allow_empty=True) if "cider_kappa" in df.columns else None
    )

    remaining_mask = pd.Series(True, index=df.index)
    filter_report = {}

    if gravy_max is not None:
        remaining_mask, report = apply_upper_bound_filter(df, remaining_mask, "gravy", gravy_max, "gravy")
        filter_report.update(report)

    if instability_max is not None:
        remaining_mask, report = apply_upper_bound_filter(
            df, remaining_mask, "instability_index", instability_max, "instability_index"
        )
        filter_report.update(report)

    if helix_sheet_max is not None:
        remaining_mask, report = apply_upper_bound_filter(
            df, remaining_mask, "helix_sheet_sum", helix_sheet_max, "helix_sheet_sum"
        )
        filter_report.update(report)

    if turn_fraction_min is not None:
        remaining_mask, report = apply_lower_bound_filter(
            df, remaining_mask, "turn_fraction", turn_fraction_min, "turn_fraction"
        )
        filter_report.update(report)

    if charge_sum_min is not None:
        remaining_mask, report = apply_lower_bound_filter(
            df, remaining_mask, "charge_sum", charge_sum_min, "charge_sum"
        )
        filter_report.update(report)

    if charge_pH7_min is not None and charge_pH7_max is not None:
        remaining_mask, report = apply_range_filter(
            df, remaining_mask, "charge_at_pH7", charge_pH7_min, charge_pH7_max, "charge_at_pH7"
        )
        filter_report.update(report)

    if fcr_max is not None:
        remaining_mask, report = apply_upper_bound_filter(df, remaining_mask, "fcr", fcr_max, "fcr")
        filter_report.update(report)

    if ncpr_min is not None and ncpr_max is not None:
        remaining_mask, report = apply_range_filter(
            df, remaining_mask, "ncpr", ncpr_min, ncpr_max, "ncpr"
        )
        filter_report.update(report)

    if cider_kappa_max is not None:
        remaining_mask, report = apply_upper_bound_filter(
            df, remaining_mask, "cider_kappa", cider_kappa_max, "cider_kappa"
        )
        filter_report.update(report)

    filtered_df = df.loc[remaining_mask].copy()
    rows_removed = len(df) - len(filtered_df)
    print(f"\nFiltered rows: {len(filtered_df):,} / {len(df):,}")
    print(f"Removed rows: {rows_removed:,} of {len(df):,}")

    if filter_report:
        print("\nPer-parameter removal:")
        for name, report in filter_report.items():
            print(
                f"  - {name}: removed {report['removed']:,}; "
                f"remaining {report['remaining_after']:,}"
            )

    output_file = output_dir / f"filtered_{input_file.stem}.csv"
    metadata = {
        "source_file": str(input_file.name),
        "rows_before": len(df),
        "rows_after": len(filtered_df),
        "rows_removed": rows_removed,
        "missing_optional_columns": ", ".join(missing_optional) if missing_optional else "",
        "skipped_filters": ", ".join(skipped_filters) if skipped_filters else "",
        "gravy_max": gravy_max,
        "instability_index_max": instability_max,
        "helix_sheet_sum_max": helix_sheet_max,
        "turn_fraction_min": turn_fraction_min,
        "charge_sum_min": charge_sum_min,
        "charge_at_pH7_min": charge_pH7_min,
        "charge_at_pH7_max": charge_pH7_max,
        "fcr_max": fcr_max,
        "ncpr_min": ncpr_min,
        "ncpr_max": ncpr_max,
        "cider_kappa_max": cider_kappa_max,
    }
    for key, report in filter_report.items():
        metadata[f"removed_by_{key}"] = report["removed"]
        metadata[f"remaining_after_{key}"] = report["remaining_after"]

    write_output_with_metadata(output_file, filtered_df, metadata)

    print(f"Filtered results written to: {output_file}")
    print("Filter parameters saved inside output CSV as metadata comments.")
    print(f"Length summary stats file: {summary_csv_path}")


if __name__ == "__main__":
    main()
