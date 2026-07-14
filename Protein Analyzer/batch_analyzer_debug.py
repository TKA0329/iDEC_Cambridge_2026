#!/usr/bin/env python3
"""
Debug version to test if progress callback is being called
"""

import sys
import time
import os
import pandas as pd
from protein_pipeline import analyze_sequences_parallel


if __name__ == "__main__":
    input_file = "sequences_only.csv"
    output_file = "results_debug.csv"

    print(f"📖 Reading {input_file}...")
    df = pd.read_csv(input_file, comment="#")

    seq_col = df.columns[0]
    print(f"✓ Found {len(df):,} sequences")

    num_cores = max(1, os.cpu_count() - 1)
    print(f"🧬 Using {num_cores} cores...")
    print(f"💡 Callback test starting...\n")

    start_time = time.time()
    callback_count = 0

    def progress_callback(completed, total):
        global callback_count
        callback_count += 1
        elapsed = time.time() - start_time
        print(f"[DEBUG] Callback #{callback_count}: {completed:,}/{total:,} | {elapsed:.1f}s elapsed")
        sys.stdout.flush()

    print("Starting analysis...")
    results = analyze_sequences_parallel(
        df[seq_col].astype(str).tolist(),
        progress_callback=progress_callback
    )

    elapsed = time.time() - start_time
    print(f"\n✓ Done! {callback_count} callbacks received in {elapsed:.1f}s")
