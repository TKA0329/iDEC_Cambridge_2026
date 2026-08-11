#!/usr/bin/env python3
"""
Batch Protein Property Analyzer — CLI Version
Run: python batch_analyzer.py input.csv output.csv
"""

import argparse
import sys
import time
import os

def main():
    # Track initialization phases
    init_start = time.time()
    print("[init] Initializing...", end="", flush=True)

    import pandas as pd
    from protein_pipeline import (
        HAS_LOCALCIDER,
        read_csv_safe,
        analyze_sequences_parallel,
        ensure_analysis_columns,
        merge_analysis_results,
        screen_analysis_columns,
    )

    init_imports = time.time() - init_start
    print(f" imports {init_imports:.2f}s", flush=True)

    try:
        import psutil
        HAS_PSUTIL = True
    except ImportError:
        HAS_PSUTIL = False
        psutil = None

    parser = argparse.ArgumentParser(description="Batch protein property analyzer")
    parser.add_argument("input_file", help="Input CSV file")
    parser.add_argument("output_file", help="Output CSV file")
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=None,
        help="Worker count. Defaults to 8 threads or CPU count for processes.",
    )
    parser.add_argument("--processes", action="store_true",
                        help="Use multiprocessing instead of threading (slower on Windows, try if you have GIL issues)")
    parser.add_argument(
        "--cider-kappa",
        action="store_true",
        help="Include localCIDER kappa in the output (slowest metric; disabled by default).",
    )
    parser.add_argument(
        "--no-fcr-ncpr",
        action="store_true",
        help="Skip FCR and NCPR calculations.",
    )
    parser.add_argument(
        "--cahs-motif1",
        action="store_true",
        help=(
            "Score the CAHS motif-1 helical segment: Pace & Scholtz helix "
            "propensity (sum/mean), Eisenberg hydrophobic moment, helical "
            "face occupancy, salt-bridge count, and Pro/Gly helix-breaker flags."
        ),
    )
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    use_processes = args.processes
    include_cider_kappa = args.cider_kappa
    include_extended_charge_metrics = not args.no_fcr_ncpr
    include_cahs_motif1 = args.cahs_motif1

    default_workers = (os.cpu_count() or 4) if use_processes else min(os.cpu_count() or 4, 8)
    num_workers = max(1, args.workers or default_workers)

    print(f"[read] Reading {input_file}...", end="", flush=True)
    read_start = time.time()
    try:
        df = read_csv_safe(input_file, comment="#")
    except FileNotFoundError:
        print(f"\n[error] File not found: {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"\n[error] Error reading file: {e}")
        sys.exit(1)
    read_elapsed = time.time() - read_start
    print(f" {read_elapsed:.2f}s", flush=True)
    
    if df.empty:
        print("[error] Input CSV is empty")
        sys.exit(1)
    
    # Get sequence column (first column)
    seq_col = df.columns[0]
    print(f"[ok] Found {len(df):,} sequences in column '{seq_col}'")
    
    if not HAS_PSUTIL:
        print("[tip] Install psutil for CPU usage stats: pip install psutil")

    if include_cider_kappa and not HAS_LOCALCIDER:
        print("[warn] localcider is not installed, so `cider_kappa` will be blank for all rows.")

    metrics_summary = ["core ProtParam metrics", "net charge"]
    if include_extended_charge_metrics:
        metrics_summary.append("FCR/NCPR")
    if include_cider_kappa:
        metrics_summary.append("CIDER kappa")
    if include_cahs_motif1:
        metrics_summary.append("CAHS motif-1 helix/amphipathicity")
    print(f"[metrics] Enabled: {', '.join(metrics_summary)}")
    
    analysis_plan = screen_analysis_columns(
        df,
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )
    rows_to_update_mask = analysis_plan["rows_to_update_mask"]
    rows_to_update_count = analysis_plan["rows_to_update_count"]

    if analysis_plan["missing_columns"]:
        print(f"[screen] Missing analyzer columns: {', '.join(analysis_plan['missing_columns'])}")
    if analysis_plan["partially_empty_columns"]:
        print(f"[screen] Columns with blank cells: {', '.join(analysis_plan['partially_empty_columns'])}")

    if rows_to_update_count == 0:
        print("[screen] All requested analyzer columns are already populated. Nothing to update.")
        out_df = ensure_analysis_columns(
            df,
            include_extended_charge_metrics=include_extended_charge_metrics,
            include_cider_kappa=include_cider_kappa,
            include_cahs_motif1=include_cahs_motif1,
        )
        results_df = out_df[analysis_plan["requested_columns"]].copy()
        elapsed = 0.0
        print(f"[write] Writing results to {output_file}...")
        try:
            out_df.to_csv(output_file, index=False)
        except Exception as e:
            print(f"[error] Error writing output file: {e}")
            sys.exit(1)

        n_ok = (results_df["error"].fillna("") == "").sum()
        n_err = len(results_df) - n_ok
        print()
        print("=" * 70)
        print("[ok] Analysis complete!")
        print(f"  Sequences in file: {len(df):,}")
        print(f"  Rows updated: 0")
        print(f"  Successful: {n_ok:,}")
        print(f"  Errors: {n_err:,}")
        print(f"  Output: {output_file}")
        print("=" * 70)
        return

    # Analyze sequences with progress tracking
    executor_type = "processes" if use_processes else "threads"
    print(
        f"[run] Updating {rows_to_update_count:,} of {len(df):,} sequences "
        f"using {num_workers} {executor_type}..."
    )
    if include_cider_kappa and HAS_LOCALCIDER:
        print("[run] Starting workers... (first few sequences may be slower due to localcider warmup)")
    else:
        print("[run] Starting workers...")
    sys.stdout.flush()
    
    # Convert sequences to list
    convert_start = time.time()
    sequences = df.loc[rows_to_update_mask, seq_col].astype(str)
    convert_elapsed = time.time() - convert_start
    
    total_init_elapsed = time.time() - init_start
    print(f"[init] {total_init_elapsed:.2f}s (imports {init_imports:.2f}s, read {read_elapsed:.2f}s, convert {convert_elapsed:.2f}s)")
    
    # Print initial 0% bar immediately
    bar_length = 35
    bar = "." * bar_length
    print(f"[{bar}]   0.00% | 0/{len(sequences):,} | CPU:   0.0% | MEM:    0MB |    0.0s |    0.0s rem")
    sys.stdout.flush()
    
    start_time = time.time()
    last_update = start_time
    process = psutil.Process() if HAS_PSUTIL else None
    update_interval = 1.0 if include_cider_kappa else 0.2

    # Prime psutil so later non-blocking reads are meaningful.
    if HAS_PSUTIL:
        try:
            psutil.cpu_percent(interval=None)
        except Exception:
            pass
    
    def progress_callback(completed, total):
        nonlocal last_update
        now = time.time()
        # Show updates periodically, always show first completion and final.
        if now - last_update < update_interval and completed != 1 and completed != total:
            return
        last_update = now
        
        elapsed = now - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        remaining = (total - completed) / rate if rate > 0 else 0
        percent = (completed / total) * 100
        bar_length = 35
        filled = int(bar_length * completed / total)
        bar = "#" * filled + "." * (bar_length - filled)
        
        # CPU usage
        cpu_str = ""
        mem_str = ""
        if HAS_PSUTIL:
            try:
                cpu_percent = psutil.cpu_percent(interval=None)
                mem_info = process.memory_info()
                mem_mb = mem_info.rss / (1024 * 1024)
                cpu_str = f" | CPU: {cpu_percent:5.1f}%"
                mem_str = f" | MEM: {mem_mb:6.0f}MB"
            except:
                pass
        
        msg = f"\r[{bar}] {percent:6.2f}% | {completed:,}/{total:,}{cpu_str}{mem_str} | {elapsed:6.1f}s | {remaining:6.1f}s rem"
        sys.stdout.write(msg)
        sys.stdout.flush()
    
    results = analyze_sequences_parallel(
        sequences,
        num_workers=num_workers,
        progress_callback=progress_callback,
        use_processes=use_processes,
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )
    
    elapsed = time.time() - start_time
    print(f"\n[ok] Analysis completed in {elapsed:.1f}s")
    
    # Merge results back into existing analyzer columns without duplicating them.
    results_df = pd.DataFrame(results, index=df.index[rows_to_update_mask])
    out_df = merge_analysis_results(
        df,
        results_df,
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )
    out_df = ensure_analysis_columns(
        out_df,
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )
    results_df = out_df[analysis_plan["requested_columns"]].copy()
    
    # Save results
    print(f"[write] Writing results to {output_file}...")
    try:
        out_df.to_csv(output_file, index=False)
    except Exception as e:
        print(f"[error] Error writing output file: {e}")
        sys.exit(1)
    
    # Summary stats
    n_ok = results_df["error"].fillna("").eq("").sum()
    n_err = len(results_df) - n_ok
    
    print()
    print("=" * 70)
    print("[ok] Analysis complete!")
    print(f"  Sequences in file: {len(df):,}")
    print(f"  Rows updated: {rows_to_update_count:,}")
    print(f"  Successful: {n_ok:,}")
    print(f"  Errors: {n_err:,}")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print(f"  Speed: {rows_to_update_count/elapsed:.0f} sequences/sec")
    print(f"  Output: {output_file}")
    print("=" * 70)


if __name__ == "__main__":
    main()