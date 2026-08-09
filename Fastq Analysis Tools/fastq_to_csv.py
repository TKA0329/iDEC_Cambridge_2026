#!/usr/bin/env python3
"""
fastq_to_csv.py

Converts a FASTQ file into a CSV of per-read statistics.

FASTQ format (4 lines per read):
    Line 1: @read_id                  (header, starts with '@')
    Line 2: SEQUENCE                  (A/T/C/G/N bases)
    Line 3: +                         (separator, may repeat the header)
    Line 4: QUALITY STRING            (same length as sequence, one char per base)

Quality encoding: Phred+33 (standard for Nanopore/Illumina 1.8+, which is what
Plasmidsaurus uses). Score for a character c is: ord(c) - 33.

Usage:
    python fastq_to_csv.py input.fastq output.csv

Requires: tqdm (pip install tqdm)
"""

import argparse
import csv
import sys

from tqdm import tqdm


def phred_scores(quality_string):
    """Convert a Phred+33 quality string into a list of integer Q scores."""
    return [ord(ch) - 33 for ch in quality_string]


def longest_homopolymer(seq):
    """Length of the longest run of a single repeated base (e.g. 'AAAA' -> 4)."""
    if not seq:
        return 0
    longest = 1
    current = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def count_reads(filepath):
    """Quickly count the number of reads in a FASTQ file (total lines / 4)."""
    with open(filepath, "r") as f:
        line_count = sum(1 for _ in f)
    return line_count // 4


def parse_fastq(filepath):
    """
    Generator that yields (read_id, sequence, quality_string) tuples,
    one per read, reading the file 4 lines at a time.
    """
    with open(filepath, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break  # end of file
            seq = f.readline().rstrip("\n")
            plus = f.readline()
            qual = f.readline().rstrip("\n")

            header = header.rstrip("\n")
            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ: expected '@' header, got: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ: expected '+' separator, got: {plus!r}")
            if len(seq) != len(qual):
                raise ValueError(
                    f"Malformed FASTQ: sequence/quality length mismatch for read {header!r}"
                )

            read_id = header[1:]  # strip leading '@'
            yield read_id, seq, qual


def compute_stats(read_id, seq, qual):
    """Compute the stats dict for a single read."""
    length = len(seq)

    a = seq.count("A")
    t = seq.count("T")
    c = seq.count("C")
    g = seq.count("G")
    n = seq.count("N")

    def pct(count):
        return round(100 * count / length, 2) if length else 0.0

    gc_pct = round(100 * (g + c) / length, 2) if length else 0.0

    scores = phred_scores(qual)
    mean_q = round(sum(scores) / len(scores), 2) if scores else 0.0
    min_q = min(scores) if scores else 0
    max_q = max(scores) if scores else 0

    # Estimated probability of error per base, derived from mean Q:
    # P(error) = 10^(-Q/10). Expressed as a percent for readability.
    est_error_rate_pct = round(100 * (10 ** (-mean_q / 10)), 4) if scores else 0.0

    return {
        "read_id": read_id,
        "seq": seq,
        "quality": qual,
        "length": length,
        "A_count": a,
        "T_count": t,
        "C_count": c,
        "G_count": g,
        "N_count": n,
        "A_pct": pct(a),
        "T_pct": pct(t),
        "C_pct": pct(c),
        "G_pct": pct(g),
        "GC_pct": gc_pct,
        "mean_quality": mean_q,
        "min_quality": min_q,
        "max_quality": max_q,
        "est_error_rate_pct": est_error_rate_pct,
        "longest_homopolymer": longest_homopolymer(seq),
    }


FIELDNAMES = [
    "read_id", "seq", "quality", "length",
    "A_count", "T_count", "C_count", "G_count", "N_count",
    "A_pct", "T_pct", "C_pct", "G_pct", "GC_pct",
    "mean_quality", "min_quality", "max_quality",
    "est_error_rate_pct", "longest_homopolymer",
]


def main():
    parser = argparse.ArgumentParser(
        description="Convert a FASTQ file into a per-read stats CSV."
    )
    parser.add_argument("input_fastq", help="Path to the input .fastq file")
    parser.add_argument("output_csv", help="Path to write the output .csv file")
    args = parser.parse_args()

    read_count = 0
    try:
        total_reads = count_reads(args.input_fastq)
    except FileNotFoundError:
        print(f"Error: could not find input file '{args.input_fastq}'", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.output_csv, "w", newline="") as out_f:
            writer = csv.DictWriter(out_f, fieldnames=FIELDNAMES)
            writer.writeheader()
            for read_id, seq, qual in tqdm(
                parse_fastq(args.input_fastq),
                total=total_reads,
                desc="Processing reads",
                unit="read",
            ):
                writer.writerow(compute_stats(read_id, seq, qual))
                read_count += 1
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Done. Wrote stats for {read_count} read(s) to '{args.output_csv}'.")


if __name__ == "__main__":
    main()