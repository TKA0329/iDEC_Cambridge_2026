#!/usr/bin/env python3
"""
CSV Chunker
===========
Splits a single-column sequences CSV into chunks of a specified max size.
Each chunk retains the 'sequence' header row.

Output files are named:
    sequences_chunk_001.csv
    sequences_chunk_002.csv
    ...

Usage:
    python3 chunk_sequences.py
    python3 chunk_sequences.py input.csv
    python3 chunk_sequences.py input.csv 200
    python3 chunk_sequences.py input.csv 200 output_folder/
"""

import csv
import os
import sys
import time

# ================================================================
#  CONFIGURATION
# ================================================================

INPUT_FILE   = 'sequences_only.csv'    # change or pass as argument
MAX_MB       = 150                   # maximum chunk size in MB
OUTPUT_DIR   = '.'                     # folder to write chunks into
CHUNK_PREFIX = 'sequences_chunk'       # output filename prefix

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

# ================================================================
#  HELPERS
# ================================================================

def _count_rows(filepath: str) -> int:
    print("  Counting rows …", end='', flush=True)
    with open(filepath, 'r', newline='') as f:
        count = sum(1 for _ in f) - 1   # subtract header
    print(f" {count:,}")
    return count


def _chunk_path(output_dir: str, index: int) -> str:
    return os.path.join(output_dir, f"{CHUNK_PREFIX}_{index:03d}.csv")


# ================================================================
#  MAIN
# ================================================================

def chunk(input_file: str, max_mb: float, output_dir: str) -> None:
    sep      = '=' * 60
    max_bytes = int(max_mb * 1024 * 1024)

    print(f"\n{sep}")
    print("  CSV Chunker")
    print(sep)
    print(f"\n  Input      : {input_file}")
    print(f"  Max size   : {max_mb} MB  ({max_bytes:,} bytes)")
    print(f"  Output dir : {output_dir}\n")

    if not os.path.exists(input_file):
        print(f"  ✗  File not found: {input_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    total = _count_rows(input_file)
    print()

    header        = ['sequence']
    header_bytes  = len(('sequence\n').encode('utf-8'))   # size of header row

    chunk_index   = 1
    chunk_bytes   = header_bytes
    chunks_done   = 0
    written_total = 0
    rows_in_chunk = 0
    t0            = time.time()

    out_path = _chunk_path(output_dir, chunk_index)
    fout     = open(out_path, 'w', newline='')
    writer   = csv.writer(fout)
    writer.writerow(header)

    with open(input_file, 'r', newline='') as fin:
        reader = csv.reader(fin)

        # Validate header
        file_header = next(reader, None)
        if file_header is None or 'sequence' not in file_header:
            print(f"  ✗  No 'sequence' column found. Got: {file_header}")
            fout.close()
            sys.exit(1)

        seq_col = file_header.index('sequence')

        for row in reader:
            seq       = row[seq_col]
            row_bytes = len((seq + '\n').encode('utf-8'))

            # If adding this row would exceed the limit, close and open next chunk
            if chunk_bytes + row_bytes > max_bytes and rows_in_chunk > 0:
                fout.close()
                chunks_done   += 1
                chunk_index   += 1
                chunk_bytes    = header_bytes
                rows_in_chunk  = 0
                out_path       = _chunk_path(output_dir, chunk_index)
                fout           = open(out_path, 'w', newline='')
                writer         = csv.writer(fout)
                writer.writerow(header)

            writer.writerow([seq])
            chunk_bytes   += row_bytes
            rows_in_chunk += 1
            written_total += 1

            if written_total % 50_000 == 0:
                _progress(written_total, total, t0)

    fout.close()
    chunks_done += 1   # close final chunk

    _progress(written_total, total, t0)
    print()            # newline after bar

    # ── Summary ──────────────────────────────────────────────
    elapsed = time.time() - t0

    print(f"\n{sep}")
    print("  COMPLETE")
    print(sep)
    print(f"  Total sequences  : {written_total:>12,}")
    print(f"  Chunks created   : {chunks_done:>12,}")
    print(f"  Runtime          : {elapsed:>11.1f} s")
    print(f"\n  Chunk breakdown:")
    total_size = 0
    for i in range(1, chunks_done + 1):
        path  = _chunk_path(output_dir, i)
        fsize = os.path.getsize(path)
        total_size += fsize
        rows  = sum(1 for _ in open(path)) - 1
        print(f"    {os.path.basename(path)}  —  "
              f"{rows:>10,} sequences  {fsize/1e6:>7.2f} MB")
    print(f"\n  Total output size : {total_size/1e6:.1f} MB")
    print(sep + "\n")


if __name__ == '__main__':
    args  = sys.argv[1:]
    inp   = args[0] if len(args) > 0 else INPUT_FILE
    mb    = float(args[1]) if len(args) > 1 else MAX_MB
    odir  = args[2] if len(args) > 2 else OUTPUT_DIR
    chunk(inp, mb, odir)
