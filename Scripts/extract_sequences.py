#!/usr/bin/env python3
"""
Sequence Extractor
==================
Reads the full provenance CSV produced by idp_generator.py and writes
a clean single-column CSV containing only the resolved sequences.

Output format:
    sequence
    MKSQTEDGP...
    ESNQTDKRP...
    ...

Usage:
    python3 extract_sequences.py
    python3 extract_sequences.py input.csv output.csv
"""

import csv
import os
import sys
import time

# ================================================================
#  CONFIGURATION
# ================================================================

INPUT_FILE  = 'generated_sequences.csv'   # change or pass as argument
OUTPUT_FILE = 'sequences_only.csv'        # change or pass as argument

# ================================================================
#  PROGRESS BAR
# ================================================================

def _progress(current: int, total: int, t0: float, bar_width: int = 40) -> None:
    elapsed = time.time() - t0
    rate    = current / elapsed if elapsed > 0 else 0
    pct     = min(current / total, 1.0) if total > 0 else 0
    eta     = (total - current) / rate if rate > 0 and total > 0 else 0
    fill    = int(bar_width * pct)
    bar     = '█' * fill + '░' * (bar_width - fill)
    print(
        f"\r  [{bar}] {pct*100:>5.1f}%  "
        f"{current:>10,}/{total:,}  "
        f"{rate:>8,.0f} rows/s  "
        f"ETA {eta:>5.0f}s",
        end='', flush=True,
    )


def _count_rows(filepath: str) -> int:
    """Count data rows (excluding header) efficiently."""
    print("  Counting rows …", end='', flush=True)
    with open(filepath, 'r', newline='') as f:
        count = sum(1 for _ in f) - 1   # subtract header
    print(f" {count:,}")
    return count


# ================================================================
#  MAIN
# ================================================================

def extract(input_file: str, output_file: str) -> None:
    sep = '=' * 60
    print(f"\n{sep}")
    print("  Sequence Extractor")
    print(sep)
    print(f"\n  Input  : {input_file}")
    print(f"  Output : {output_file}\n")

    if not os.path.exists(input_file):
        print(f"  ✗  File not found: {input_file}")
        sys.exit(1)

    total = _count_rows(input_file)
    print()

    written = 0
    t0      = time.time()

    with open(input_file,  'r', newline='') as fin, \
         open(output_file, 'w', newline='') as fout:

        reader = csv.DictReader(fin)

        if 'sequence' not in (reader.fieldnames or []):
            print(f"  ✗  No 'sequence' column found in {input_file}.")
            print(f"     Available columns: {reader.fieldnames}")
            sys.exit(1)

        writer = csv.writer(fout)
        writer.writerow(['sequence'])       # header

        for row in reader:
            writer.writerow([row['sequence']])
            written += 1
            if written % 50_000 == 0:
                _progress(written, total, t0)

    _progress(written, total, t0)          # final bar at 100 %
    print()                                # newline after bar

    elapsed = time.time() - t0
    fsize   = os.path.getsize(output_file)

    print(f"\n{sep}")
    print("  COMPLETE")
    print(sep)
    print(f"  Sequences written : {written:>12,}")
    print(f"  Runtime           : {elapsed:>11.1f} s")
    print(f"  Output file       : {output_file}")
    print(f"  File size         : {fsize/1e6:>11.1f} MB")
    print(sep + "\n")


if __name__ == '__main__':
    args = sys.argv[1:]
    inp  = args[0] if len(args) > 0 else INPUT_FILE
    out  = args[1] if len(args) > 1 else OUTPUT_FILE
    extract(inp, out)
