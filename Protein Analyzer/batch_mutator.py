#!/usr/bin/env python3
"""
Batch Protein Mutator — CLI Version

Generates mutant sequences for every row in an input CSV, using one of three
mutation modes:

  random       — each position in the region is replaced by a random AA
                 (excluding the original)
  conservative — each position is replaced by a biochemically-similar AA
                 drawn from the CONSERVATIVE_GROUPS table
  scan         — exhaustive single-position scan: one mutant per conservative
                 substitute per position (copies_col is ignored)

Run:
    python batch_mutator.py input.csv output.csv
    python batch_mutator.py input.csv output.csv --mode conservative -j 8
    python batch_mutator.py input.csv output.csv --mode scan --seq-col sequence
"""

import argparse
import random
import sys
import time
import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

import pandas as pd

# ── Import from your existing pipeline ──────────────────────────────────────
try:
    from protein_pipeline import (
        CONSERVATIVE_GROUPS,
        VALID_AA,
        parse_mutation_regions,
        parse_copy_count,
        max_unique_conservative_variants,
        count_random_combos,
        scan_sequence_single_position,
        mutate_sequence_random,
        mutate_sequence_conservative,
    )
except ImportError as e:
    print(f"[error] Could not import protein_pipeline: {e}")
    print("        Make sure protein_pipeline.py is in the same directory.")
    sys.exit(1)

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# ── Column auto-detection ─────────────────────────────────────────────────────

_REGION_KEYWORDS  = ["region", "mut_region", "mutation_region", "mutregion", "positions"]
_COPIES_KEYWORDS  = ["copies", "copy", "num_copies", "n_copies", "count", "replicates","num_random_copies"]


def _find_column(df_columns, candidates, label):
    lower_cols = {c.lower(): c for c in df_columns}
    for kw in candidates:
        if kw in lower_cols:
            return lower_cols[kw], "exact"
    for kw in candidates:
        for lc, orig in lower_cols.items():
            if kw in lc:
                return orig, "fuzzy"
    return None, None


# ── Per-row worker functions ──────────────────────────────────────────────────
#
# PROGRESS TRACKING DESIGN
# ────────────────────────
# Each worker receives a shared threading.Value (or multiprocessing.Value)
# counter called `progress_counter`.  Workers increment it by 1 for every
# mutation copy they *finish generating*, not just once per row.  The main
# thread polls that counter on a short interval to render a smooth bar that
# reflects work done at the copy level, not the row level.
#
# Total units = sum of num_copies across all rows (pre-computed in main()).
# For scan mode, total units = total number of variants generated.
#
# Workers also receive `row_total` — the number of copies this row will emit —
# so the main thread can show a "current row" sub-count if desired.
#
# NOTE: When using multiprocessing (--processes), Python's
# multiprocessing.Value is used instead of threading.Value; both expose the
# same .value attribute and a context-manager lock, so the workers below are
# identical for both executors.

def _process_row_random(args):
    """Process a single row in random mutation mode."""
    idx, row_dict, seq_col, region_col, copies_col, seed, progress_counter = args
    row = pd.Series(row_dict)
    raw_seq = row[seq_col]
    raw_region = row.get(region_col, "")
    raw_copies = row.get(copies_col, 1)
    region_text = "" if pd.isna(raw_region) else str(raw_region).strip()

    # ── Pass through row unchanged if region or copies cell is blank ───────────
    if not region_text or pd.isna(raw_copies) or str(raw_copies).strip() == "":
        passthrough = dict(row_dict)
        passthrough["source_row"] = idx + 1
        passthrough["mutation_mode"] = "skipped"
        passthrough["mutation_error"] = ""
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += 1
        return [passthrough], True, None

    rng = random.Random(seed)
    row_error = ""

    try:
        regions = parse_mutation_regions(region_text)
    except ValueError as e:
        regions = []
        row_error = str(e)

    try:
        num_copies = parse_copy_count(raw_copies)
    except ValueError as e:
        num_copies = 1
        row_error = str(e) if not row_error else f"{row_error}; {e}"

    if num_copies == 0:
        return [], True, None

    dup_warning = None
    if not row_error and regions:
        try:
            max_unique, n_positions = count_random_combos(str(raw_seq), regions)
            if num_copies > max_unique:
                dup_warning = {
                    "row": idx + 1,
                    "requested": num_copies,
                    "max_unique": max_unique,
                    "mode": "random",
                }
        except Exception:
            pass

    expanded = []
    for copy_i in range(1, num_copies + 1):
        new_row = dict(row_dict)
        new_row["source_row"] = idx + 1
        new_row["copy_index"] = copy_i
        new_row["requested_copies"] = num_copies
        new_row["original_sequence"] = str(raw_seq)
        new_row["mutation_regions"] = region_text
        new_row["mutation_mode"] = "random"

        if row_error:
            new_row[seq_col] = str(raw_seq)
            new_row["mutated_positions_count"] = 0
            new_row["conservative_skipped"] = ""
            new_row["mutation_error"] = row_error
        else:
            try:
                clean = str(raw_seq).strip().upper().replace(" ", "").replace("\n", "")
                bad = set(clean) - set(VALID_AA)
                if bad:
                    raise ValueError(f"invalid chars: {bad}")
                if not clean:
                    raise ValueError("empty sequence")

                arr = list(clean)
                n_changed = 0
                for start, end in regions:
                    for pos in range(start, end + 1):
                        i = pos - 1
                        if 0 <= i < len(arr):
                            choices = [aa for aa in VALID_AA if aa != arr[i]]
                            arr[i] = rng.choice(choices)
                            n_changed += 1

                new_row[seq_col] = "".join(arr)
                new_row["mutated_positions_count"] = n_changed
                new_row["conservative_skipped"] = ""
                new_row["mutation_error"] = ""
            except ValueError as e:
                new_row[seq_col] = str(raw_seq)
                new_row["mutated_positions_count"] = 0
                new_row["conservative_skipped"] = ""
                new_row["mutation_error"] = str(e)

        expanded.append(new_row)

        # ── Intra-row progress: increment after each copy is done ──────────
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += 1

    return expanded, False, dup_warning


def _process_row_conservative(args):
    """Process a single row in conservative mutation mode."""
    idx, row_dict, seq_col, region_col, copies_col, seed, progress_counter = args
    row = pd.Series(row_dict)
    raw_seq = row[seq_col]
    raw_region = row.get(region_col, "")
    raw_copies = row.get(copies_col, 1)
    region_text = "" if pd.isna(raw_region) else str(raw_region).strip()

    # ── Pass through row unchanged if region or copies cell is blank ───────────
    if not region_text or pd.isna(raw_copies) or str(raw_copies).strip() == "":
        passthrough = dict(row_dict)
        passthrough["source_row"] = idx + 1
        passthrough["mutation_mode"] = "skipped"
        passthrough["mutation_error"] = ""
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += 1
        return [passthrough], True, None

    rng = random.Random(seed)
    row_error = ""

    try:
        regions = parse_mutation_regions(region_text)
    except ValueError as e:
        regions = []
        row_error = str(e)

    try:
        num_copies = parse_copy_count(raw_copies)
    except ValueError as e:
        num_copies = 1
        row_error = str(e) if not row_error else f"{row_error}; {e}"

    if num_copies == 0:
        return [], True, None

    dup_warning = None
    if not row_error and regions:
        try:
            max_unique, per_pos = max_unique_conservative_variants(str(raw_seq), regions)
            if num_copies > max_unique:
                dup_warning = {
                    "row": idx + 1,
                    "requested": num_copies,
                    "max_unique": max_unique,
                    "mode": "conservative",
                    "per_position": per_pos,
                }
        except Exception:
            pass

    expanded = []
    for copy_i in range(1, num_copies + 1):
        new_row = dict(row_dict)
        new_row["source_row"] = idx + 1
        new_row["copy_index"] = copy_i
        new_row["requested_copies"] = num_copies
        new_row["original_sequence"] = str(raw_seq)
        new_row["mutation_regions"] = region_text
        new_row["mutation_mode"] = "conservative"

        if row_error:
            new_row[seq_col] = str(raw_seq)
            new_row["mutated_positions_count"] = 0
            new_row["conservative_skipped"] = ""
            new_row["mutation_error"] = row_error
        else:
            try:
                clean = str(raw_seq).strip().upper().replace(" ", "").replace("\n", "")
                bad = set(clean) - set(VALID_AA)
                if bad:
                    raise ValueError(f"invalid chars: {bad}")
                if not clean:
                    raise ValueError("empty sequence")

                arr = list(clean)
                n_changed = 0
                skipped = []
                for start, end in regions:
                    for pos in range(start, end + 1):
                        i = pos - 1
                        if 0 <= i < len(arr):
                            choices = CONSERVATIVE_GROUPS.get(arr[i], [])
                            if choices:
                                arr[i] = rng.choice(choices)
                                n_changed += 1
                            else:
                                skipped.append(pos)

                new_row[seq_col] = "".join(arr)
                new_row["mutated_positions_count"] = n_changed
                new_row["conservative_skipped"] = ",".join(map(str, skipped)) if skipped else ""
                new_row["mutation_error"] = ""
            except ValueError as e:
                new_row[seq_col] = str(raw_seq)
                new_row["mutated_positions_count"] = 0
                new_row["conservative_skipped"] = ""
                new_row["mutation_error"] = str(e)

        expanded.append(new_row)

        # ── Intra-row progress: increment after each copy is done ──────────
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += 1

    return expanded, False, dup_warning


def _process_row_scan(args):
    """Process a single row in exhaustive scan mode."""
    idx, row_dict, seq_col, region_col, progress_counter = args
    row = pd.Series(row_dict)
    raw_seq = row[seq_col]
    raw_region = row.get(region_col, "")
    region_text = "" if pd.isna(raw_region) else str(raw_region).strip()

    # ── Pass through row unchanged if region cell is blank ─────────────────────
    if not region_text:
        passthrough = dict(row_dict)
        passthrough["source_row"] = idx + 1
        passthrough["mutation_mode"] = "skipped"
        passthrough["mutation_error"] = ""
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += 1
        return [passthrough], True, None

    base_row = dict(row_dict)
    base_row["source_row"] = idx + 1
    base_row["original_sequence"] = str(raw_seq)
    base_row["mutation_regions"] = region_text
    base_row["mutation_mode"] = "scan"

    def _bump(n=1):
        if progress_counter is not None:
            with progress_counter.get_lock():
                progress_counter.value += n

    try:
        regions = parse_mutation_regions(region_text)
    except ValueError as e:
        base_row["mutation_error"] = str(e)
        base_row["mutated_position"] = ""
        base_row["original_aa"] = ""
        base_row["new_aa"] = ""
        _bump(1)
        return [base_row], False, None

    try:
        variants = scan_sequence_single_position(raw_seq, regions)
    except ValueError as e:
        base_row["mutation_error"] = str(e)
        base_row["mutated_position"] = ""
        base_row["original_aa"] = ""
        base_row["new_aa"] = ""
        _bump(1)
        return [base_row], False, None

    if not variants:
        base_row["mutation_error"] = "no conservative substitutes found in region"
        base_row["mutated_position"] = ""
        base_row["original_aa"] = ""
        base_row["new_aa"] = ""
        _bump(1)
        return [base_row], False, None

    expanded = []
    for v in variants:
        new_row = dict(base_row)
        new_row[seq_col] = v["mutated_seq"]
        new_row["mutated_position"] = v["mutated_position"]
        new_row["original_aa"] = v["original_aa"]
        new_row["new_aa"] = v["new_aa"]
        new_row["mutation_error"] = ""
        expanded.append(new_row)
        _bump(1)  # one tick per variant generated

    return expanded, False, None


# ── Progress bar helper ──────────────────────────────────────────────────────

def _render_bar(completed, total, elapsed, process, bar_length=35):
    """
    Render a progress bar line.

    `completed` and `total` are now copy-level counts (not row counts), giving
    smooth sub-row progress when a single row generates millions of copies.
    """
    rate = completed / elapsed if elapsed > 0 else 0
    remaining = (total - completed) / rate if rate > 0 else 0
    percent = (completed / total) * 100 if total > 0 else 0
    filled = int(bar_length * completed / total) if total > 0 else 0
    bar = "#" * filled + "." * (bar_length - filled)

    cpu_str = mem_str = ""
    if HAS_PSUTIL and process is not None:
        try:
            cpu_str = f" | CPU: {psutil.cpu_percent(interval=None):5.1f}%"
            mem_str = f" | MEM: {process.memory_info().rss / (1024 * 1024):6.0f}MB"
        except Exception:
            pass

    return (
        f"\r[{bar}] {percent:6.2f}% | {completed:,}/{total:,} copies"
        f"{cpu_str}{mem_str} | {elapsed:6.1f}s | {remaining:6.1f}s rem"
    )


# ── Total-copies pre-scan ────────────────────────────────────────────────────

def _precount_total_copies(records, seq_col, region_col, copies_col, mode, default_copies):
    """
    Walk every row and sum up the number of mutation copies that will be
    generated.  This is needed so the progress bar has an accurate denominator
    at copy granularity.

    For scan mode we approximate: count conservative substitutes per position
    across the region.  This is fast (no sequence generation) and exact.

    Returns total_copies (int).  Falls back gracefully on any error.
    """
    total = 0
    for rec in records:
        raw_copies = rec.get(copies_col, default_copies)
        raw_region = rec.get(region_col, "")
        raw_seq    = rec.get(seq_col, "")
        region_text = "" if pd.isna(raw_region) else str(raw_region).strip()

        if mode == "scan":
            # Each position × number of conservative substitutes = variants
            try:
                regions = parse_mutation_regions(region_text)
                variants = scan_sequence_single_position(str(raw_seq), regions)
                total += max(1, len(variants))  # at least 1 (error row)
            except Exception:
                total += 1
        else:
            try:
                n = parse_copy_count(raw_copies)
                total += n
            except Exception:
                total += 1
    return max(total, 1)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    init_start = time.time()
    print("[init] Initializing...", flush=True)

    parser = argparse.ArgumentParser(
        description="Batch protein mutator — generates mutant sequences from a CSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("output_file", help="Output CSV file")
    parser.add_argument(
        "--mode",
        choices=["random", "conservative", "scan"],
        default="random",
        help="Mutation strategy (default: random)",
    )
    parser.add_argument(
        "--seq-col",
        default=None,
        help="Name of the sequence column (default: first column)",
    )
    parser.add_argument(
        "--region-col",
        default="mutation_regions",
        help="Name of the mutation-region column (default: mutation_regions)",
    )
    parser.add_argument(
        "--copies-col",
        default="num_copies",
        help="Name of the copy-count column, random/conservative only (default: num_copies)",
    )
    parser.add_argument(
        "--default-copies",
        type=int,
        default=1,
        help="Copies to generate when copies-col is missing/blank (default: 1)",
    )
    parser.add_argument(
        "--default-region",
        default="",
        help="Fallback region spec when region-col is missing/blank (e.g. '1-10')",
    )
    parser.add_argument(
        "-j", "--workers",
        type=int,
        default=None,
        help="Worker thread count (default: min(CPU count, 8))",
    )
    parser.add_argument(
        "--processes",
        action="store_true",
        help="Use multiprocessing instead of threading",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Base random seed (default: 42); each row gets seed+row_index",
    )
    parser.add_argument(
        "--warn-duplicates",
        action="store_true",
        help="Print a warning when requested copies exceed unique possibilities",
    )
    args = parser.parse_args()

    use_processes = args.processes
    default_workers = (os.cpu_count() or 4) if use_processes else min(os.cpu_count() or 4, 8)
    num_workers = max(1, args.workers or default_workers)
    executor_type = "processes" if use_processes else "threads"

    # ── Read input ────────────────────────────────────────────────────────────
    print(f"[read] Reading {args.input_file}...", end="", flush=True)
    read_start = time.time()
    try:
        df = pd.read_csv(
            args.input_file,
            comment="#",
            on_bad_lines="warn",   # warn and skip unparseable rows instead of crashing
            engine="python",       # python engine tolerates ragged trailing commas
        )
    except FileNotFoundError:
        print(f"\n[error] File not found: {args.input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[error] Could not read file: {e}")
        sys.exit(1)
    print(f" {time.time() - read_start:.2f}s", flush=True)

    if df.empty:
        print("[error] Input CSV is empty")
        sys.exit(1)

    # Drop phantom columns produced by trailing commas in the CSV
    # (pandas names them "Unnamed: N" when the header has fewer fields than data rows)
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    if unnamed:
        print(f"[info] Dropping {len(unnamed)} empty trailing column(s) from CSV: {unnamed}")
        df.drop(columns=unnamed, inplace=True)
    # Also drop any column that is entirely NaN (completely empty column)
    all_null = [c for c in df.columns if df[c].isna().all()]
    if all_null:
        print(f"[info] Dropping {len(all_null)} fully-empty column(s): {all_null}")
        df.drop(columns=all_null, inplace=True)

    seq_col = args.seq_col or df.columns[0]
    if seq_col not in df.columns:
        print(f"[error] Sequence column '{seq_col}' not found. Columns: {list(df.columns)}")
        sys.exit(1)

    region_col = args.region_col
    copies_col = args.copies_col
    if region_col not in df.columns:
        print(f"[info] Region column '{region_col}' not found — using default region: '{args.default_region or '(full sequence)'}'")
        df[region_col] = args.default_region
    if copies_col not in df.columns and args.mode != "scan":
        print(f"[info] Copies column '{copies_col}' not found — defaulting to {args.default_copies} copy/copies per row")
        df[copies_col] = args.default_copies

    n_rows = len(df)
    print(f"[ok] {n_rows:,} input rows | seq_col='{seq_col}' | mode={args.mode} | workers={num_workers} {executor_type}")

    if not HAS_PSUTIL:
        print("[tip] Install psutil for CPU/memory stats: pip install psutil")

    records = df.to_dict(orient="records")

    # ── Pre-count total copies for accurate progress denominator ──────────────
    print("[scan] Pre-counting total mutation copies...", end="", flush=True)
    scan_start = time.time()
    total_copies = _precount_total_copies(
        records, seq_col, region_col, copies_col, args.mode, args.default_copies
    )
    print(f" {total_copies:,} copies ({time.time() - scan_start:.2f}s)", flush=True)

    # ── Shared counter — threading.Value for threads, multiprocessing.Value for procs ──
    if use_processes:
        import multiprocessing
        progress_counter = multiprocessing.Value("q", 0)  # signed 64-bit int
    else:
        import threading
        import ctypes

        class _ThreadingValue:
            """Minimal drop-in for multiprocessing.Value using a threading lock."""
            def __init__(self, val=0):
                self.value = val
                self._lock = threading.Lock()

            def get_lock(self):
                return self._lock

        progress_counter = _ThreadingValue(0)

    # ── Build per-row task arguments ──────────────────────────────────────────
    if args.mode == "random":
        tasks = [
            (idx, rec, seq_col, region_col, copies_col, args.seed + idx, progress_counter)
            for idx, rec in enumerate(records)
        ]
        worker_fn = _process_row_random
    elif args.mode == "conservative":
        tasks = [
            (idx, rec, seq_col, region_col, copies_col, args.seed + idx, progress_counter)
            for idx, rec in enumerate(records)
        ]
        worker_fn = _process_row_conservative
    else:  # scan
        tasks = [
            (idx, rec, seq_col, region_col, progress_counter)
            for idx, rec in enumerate(records)
        ]
        worker_fn = _process_row_scan

    # ── Parallel execution with copy-level progress bar ───────────────────────
    print(f"[run] Generating mutations ({total_copies:,} total copies across {n_rows:,} rows)...")
    bar_length = 35
    print(
        f"[{'.' * bar_length}]   0.00% | 0/{total_copies:,} copies"
        + (" | CPU:   0.0% | MEM:    0MB" if HAS_PSUTIL else "")
        + " |    0.0s |    0.0s rem"
    )
    sys.stdout.flush()

    start_time = time.time()
    last_update = start_time
    update_interval = 0.2

    process = psutil.Process() if HAS_PSUTIL else None
    if HAS_PSUTIL:
        try:
            psutil.cpu_percent(interval=None)
        except Exception:
            pass

    all_expanded_rows = [None] * n_rows
    skipped_zero = 0
    dup_warnings = []
    rows_completed = 0  # tracks row-level completion for result collection

    executor_cls = ProcessPoolExecutor if use_processes else ThreadPoolExecutor

    with executor_cls(max_workers=num_workers) as executor:
        future_to_task_idx = {executor.submit(worker_fn, task): i for i, task in enumerate(tasks)}

        # ── Progress-bar polling loop ─────────────────────────────────────────
        # We use a non-blocking approach: check as_completed() with a short
        # timeout, render the bar between checks, and collect results as they
        # arrive.  This way the bar updates smoothly even while a single huge
        # row is being processed.
        pending = set(future_to_task_idx.keys())

        while pending:
            # Poll for any newly-completed futures (non-blocking snapshot)
            done_now = {f for f in pending if f.done()}

            for future in done_now:
                pending.discard(future)
                task_idx = future_to_task_idx[future]
                try:
                    expanded, was_skipped, dup_warn = future.result()
                except Exception as e:
                    rec = records[task_idx]
                    err_row = dict(rec)
                    err_row["source_row"] = task_idx + 1
                    err_row["mutation_error"] = str(e)
                    expanded = [err_row]
                    was_skipped = False
                    dup_warn = None

                all_expanded_rows[task_idx] = expanded
                if was_skipped:
                    skipped_zero += 1
                if dup_warn and args.warn_duplicates:
                    dup_warnings.append(dup_warn)
                rows_completed += 1

            # Render bar from the shared copy counter
            now = time.time()
            if now - last_update >= update_interval or not pending:
                last_update = now
                elapsed = now - start_time
                copies_done = progress_counter.value
                sys.stdout.write(_render_bar(copies_done, total_copies, elapsed, process))
                sys.stdout.flush()

            if pending:
                time.sleep(0.05)  # 50 ms sleep keeps CPU usage low between polls

    elapsed = time.time() - start_time
    # Final bar at 100 %
    sys.stdout.write(_render_bar(total_copies, total_copies, elapsed, process))
    sys.stdout.flush()
    print(f"\n[ok] All rows processed in {elapsed:.1f}s")

    # ── Duplicate warnings ────────────────────────────────────────────────────
    if dup_warnings:
        print(f"\n[warn] {len(dup_warnings)} row(s) requested more copies than unique variants:")
        for w in dup_warnings[:10]:
            print(f"       row {w['row']}: requested={w['requested']}, max_unique={w['max_unique']} ({w['mode']})")
        if len(dup_warnings) > 10:
            print(f"       ... and {len(dup_warnings) - 10} more")

    # ── Flatten & write output ────────────────────────────────────────────────
    print("[build] Assembling output dataframe...", end="", flush=True)
    flat = []
    for rows in all_expanded_rows:
        if rows:
            flat.extend(rows)
    out_df = pd.DataFrame(flat)
    print(f" {len(out_df):,} rows", flush=True)

    print(f"[write] Writing to {args.output_file}...", end="", flush=True)
    try:
        out_df.to_csv(args.output_file, index=False)
    except Exception as e:
        print(f"\n[error] Could not write output: {e}")
        sys.exit(1)
    print(" done", flush=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    n_errors = int(out_df.get("mutation_error", pd.Series(dtype=str)).fillna("").ne("").sum())
    print()
    print("=" * 70)
    print("[ok] Mutation generation complete!")
    print(f"  Input rows:         {n_rows:,}")
    print(f"  Output rows:        {len(out_df):,}")
    print(f"  Passed through (blank region/copies): {skipped_zero:,}")
    print(f"  Rows with errors:   {n_errors:,}")
    print(f"  Time elapsed:       {elapsed:.1f}s")
    print(f"  Speed:              {n_rows / elapsed:.0f} input rows/sec")
    print(f"  Output:             {args.output_file}")
    if args.warn_duplicates:
        print(f"  Duplicate warnings: {len(dup_warnings):,}")
    print("=" * 70)


if __name__ == "__main__":
    main()