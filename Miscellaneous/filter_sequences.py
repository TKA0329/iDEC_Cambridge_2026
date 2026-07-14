#!/usr/bin/env python3
"""
CSV filter for protein sequence results.

This script adds derived summary columns when possible, prints per-length stats
to the terminal, writes a length-by-length summary CSV, and generates a filtered
output CSV with filter parameters embedded as a metadata section.

Two ways to supply filter thresholds:

  1. Interactive (default when no threshold flags are given): the script
     prompts for each available filter one at a time, same as before.

  2. Non-interactive (pass --non-interactive, optionally with one or more
     of the --*-max / --*-min / --*-range flags below): no prompts are
     shown. Any filter whose flag isn't supplied is simply skipped. This
     is the mode to use when calling this script repeatedly from another
     script/pipeline with the same thresholds every time.

Run (interactive):
  python filter_sequences.py input.csv

Run (non-interactive, e.g. from another script):
  python filter_sequences.py input.csv output_dir --non-interactive \\
      --gravy-max 0.5 --instability-max 40 --ncpr-min -0.3 --ncpr-max 0.3 \\
      --pace-scholtz-sum-max 5 --pace-scholtz-mean-max 0.5 \\
      --hydrophobic-moment-min 0.2 --hydrophobic-moment-max 0.6 \\
      --helical-face-occupancy-min 0.7

Optional positional:
  python filter_sequences.py input.csv output_dir
"""

import argparse
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
    next_mask = remaining_mask & (df[column].astype(float) >= value)
    removed = int(remaining_mask.sum() - next_mask.sum())
    return next_mask, {
        report_name: {
            "value": value,
            "removed": removed,
            "remaining_after": int(next_mask.sum()),
        }
    }


def apply_lower_bound_inclusive_filter(df, remaining_mask, column, value, report_name):
    next_mask = remaining_mask & (df[column].astype(float) >= value)
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


def apply_exclude_range_filter(df, remaining_mask, column, low, high, report_name):
    """Keep rows OUTSIDE the given [low, high] range (i.e. remove rows within it)."""
    values = df[column].astype(float)
    next_mask = remaining_mask & ((values < low) | (values > high))
    removed = int(remaining_mask.sum() - next_mask.sum())
    return next_mask, {
        report_name: {
            "value": f"exclude {low}-{high}",
            "removed": removed,
            "remaining_after": int(next_mask.sum()),
        }
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a protein sequence results CSV by derived property thresholds."
    )
    parser.add_argument("input_file", type=Path, help="Input CSV file")
    parser.add_argument(
        "output_dir", type=Path, nargs="?", default=Path("filtered_output"),
        help="Directory to write outputs to (default: filtered_output)",
    )
    parser.add_argument(
        "--non-interactive", action="store_true",
        help="Do not prompt for filter values. Use the --*-max/--*-min/--*-range "
             "flags below instead; any filter whose flag is omitted is skipped.",
    )
    parser.add_argument("--gravy-max", type=float, default=None)
    parser.add_argument("--instability-max", type=float, default=None)
    parser.add_argument("--helix-sheet-max", type=float, default=None)
    parser.add_argument("--turn-fraction-min", type=float, default=None)
    parser.add_argument("--charge-sum-min", type=float, default=None)
    parser.add_argument(
        "--charge-ph7-exclude-min", type=float, default=None,
        help="Remove rows with charge_at_pH7 in [min, max]. Supply together with --charge-ph7-exclude-max.",
    )
    parser.add_argument("--charge-ph7-exclude-max", type=float, default=None)
    parser.add_argument("--fcr-max", type=float, default=None)
    parser.add_argument(
        "--ncpr-min", type=float, default=None,
        help="Keep rows with ncpr in [min, max]. Supply together with --ncpr-max.",
    )
    parser.add_argument("--ncpr-max", type=float, default=None)
    parser.add_argument("--cider-kappa-max", type=float, default=None)
    parser.add_argument("--pace-scholtz-sum-max", type=float, default=None)
    parser.add_argument("--pace-scholtz-mean-max", type=float, default=None)
    parser.add_argument(
        "--hydrophobic-moment-min", type=float, default=None,
        help="Keep rows with hydrophobic_moment in [min, max]. Supply together with --hydrophobic-moment-max.",
    )
    parser.add_argument("--hydrophobic-moment-max", type=float, default=None)
    parser.add_argument("--helical-face-occupancy-min", type=float, default=None)
    args = parser.parse_args()

    if (args.charge_ph7_exclude_min is None) != (args.charge_ph7_exclude_max is None):
        parser.error("--charge-ph7-exclude-min and --charge-ph7-exclude-max must be supplied together.")
    if (args.ncpr_min is None) != (args.ncpr_max is None):
        parser.error("--ncpr-min and --ncpr-max must be supplied together.")
    if (args.hydrophobic_moment_min is None) != (args.hydrophobic_moment_max is None):
        parser.error("--hydrophobic-moment-min and --hydrophobic-moment-max must be supplied together.")

    return args


def main():
    args = parse_args()
    input_file = args.input_file
    output_dir = args.output_dir

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
        available_filter_descriptions.append("turn_fraction: keep rows with turn_fraction >= value")
    else:
        skipped_filters.append("turn_fraction")
    if "charge_sum" in df.columns:
        available_filter_descriptions.append(
            "charge_sum: keep rows with positive_res_RK + negative_res_DE >= value"
        )
    else:
        skipped_filters.append("charge_sum")
    if "charge_at_pH7" in df.columns:
        available_filter_descriptions.append("charge_at_pH7: remove rows within the specified range")
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
    if "pace_scholtz_sum" in df.columns:
        available_filter_descriptions.append("pace_scholtz_sum: keep rows with pace_scholtz_sum <= value")
    else:
        skipped_filters.append("pace_scholtz_sum")
    if "pace_scholtz_mean" in df.columns:
        available_filter_descriptions.append("pace_scholtz_mean: keep rows with pace_scholtz_mean <= value")
    else:
        skipped_filters.append("pace_scholtz_mean")
    if "hydrophobic_moment" in df.columns:
        available_filter_descriptions.append("hydrophobic_moment: keep rows with hydrophobic_moment within range")
    else:
        skipped_filters.append("hydrophobic_moment")
    if "helical_face_occupancy" in df.columns:
        available_filter_descriptions.append(
            "helical_face_occupancy: keep rows with helical_face_occupancy >= value"
        )
    else:
        skipped_filters.append("helical_face_occupancy")

    if args.non_interactive:
        print("Non-interactive mode: using thresholds supplied on the command line "
              "(any filter without a flag is skipped).")
        gravy_max = args.gravy_max if "gravy" in df.columns else None
        instability_max = args.instability_max if "instability_index" in df.columns else None
        helix_sheet_max = args.helix_sheet_max if "helix_sheet_sum" in df.columns else None
        turn_fraction_min = args.turn_fraction_min if "turn_fraction" in df.columns else None
        charge_sum_min = args.charge_sum_min if "charge_sum" in df.columns else None
        if "charge_at_pH7" in df.columns:
            charge_pH7_min, charge_pH7_max = args.charge_ph7_exclude_min, args.charge_ph7_exclude_max
        else:
            charge_pH7_min, charge_pH7_max = None, None
        fcr_max = args.fcr_max if "fcr" in df.columns else None
        if "ncpr" in df.columns:
            ncpr_min, ncpr_max = args.ncpr_min, args.ncpr_max
        else:
            ncpr_min, ncpr_max = None, None
        cider_kappa_max = args.cider_kappa_max if "cider_kappa" in df.columns else None
        pace_scholtz_sum_max = args.pace_scholtz_sum_max if "pace_scholtz_sum" in df.columns else None
        pace_scholtz_mean_max = args.pace_scholtz_mean_max if "pace_scholtz_mean" in df.columns else None
        if "hydrophobic_moment" in df.columns:
            hydrophobic_moment_min, hydrophobic_moment_max = (
                args.hydrophobic_moment_min, args.hydrophobic_moment_max
            )
        else:
            hydrophobic_moment_min, hydrophobic_moment_max = None, None
        helical_face_occupancy_min = (
            args.helical_face_occupancy_min if "helical_face_occupancy" in df.columns else None
        )
    else:
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
            prompt_float("turn_fraction >=", allow_empty=True) if "turn_fraction" in df.columns else None
        )
        charge_sum_min = (
            prompt_float("positive_res_RK + negative_res_DE >=", allow_empty=True)
            if "charge_sum" in df.columns
            else None
        )
        if "charge_at_pH7" in df.columns:
            charge_pH7_min, charge_pH7_max = prompt_range(
                "remove charge_at_pH7 range (min-max)", allow_empty=True
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
        pace_scholtz_sum_max = (
            prompt_float("pace_scholtz_sum <=", allow_empty=True)
            if "pace_scholtz_sum" in df.columns
            else None
        )
        pace_scholtz_mean_max = (
            prompt_float("pace_scholtz_mean <=", allow_empty=True)
            if "pace_scholtz_mean" in df.columns
            else None
        )
        if "hydrophobic_moment" in df.columns:
            hydrophobic_moment_min, hydrophobic_moment_max = prompt_range(
                "hydrophobic_moment range (min-max)", allow_empty=True
            )
        else:
            hydrophobic_moment_min, hydrophobic_moment_max = None, None
        helical_face_occupancy_min = (
            prompt_float("helical_face_occupancy >=", allow_empty=True)
            if "helical_face_occupancy" in df.columns
            else None
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
        remaining_mask, report = apply_exclude_range_filter(
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

    if pace_scholtz_sum_max is not None:
        remaining_mask, report = apply_upper_bound_filter(
            df, remaining_mask, "pace_scholtz_sum", pace_scholtz_sum_max, "pace_scholtz_sum"
        )
        filter_report.update(report)

    if pace_scholtz_mean_max is not None:
        remaining_mask, report = apply_upper_bound_filter(
            df, remaining_mask, "pace_scholtz_mean", pace_scholtz_mean_max, "pace_scholtz_mean"
        )
        filter_report.update(report)

    if hydrophobic_moment_min is not None and hydrophobic_moment_max is not None:
        remaining_mask, report = apply_range_filter(
            df, remaining_mask, "hydrophobic_moment",
            hydrophobic_moment_min, hydrophobic_moment_max, "hydrophobic_moment"
        )
        filter_report.update(report)

    if helical_face_occupancy_min is not None:
        remaining_mask, report = apply_lower_bound_inclusive_filter(
            df, remaining_mask, "helical_face_occupancy", helical_face_occupancy_min, "helical_face_occupancy"
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
        "pass_rate_pct": round(100.0 * len(filtered_df) / len(df), 4) if len(df) else 0.0,
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
        "pace_scholtz_sum_max": pace_scholtz_sum_max,
        "pace_scholtz_mean_max": pace_scholtz_mean_max,
        "hydrophobic_moment_min": hydrophobic_moment_min,
        "hydrophobic_moment_max": hydrophobic_moment_max,
        "helical_face_occupancy_min": helical_face_occupancy_min,
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