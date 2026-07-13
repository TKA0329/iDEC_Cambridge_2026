#!/usr/bin/env python3
"""
filter_sequences.py

Builds a degenerate-motif representation from a set of aligned, mutated
protein sequences, based on per-position amino acid frequency.

ALGORITHM (per aligned position):
  1. Count each amino acid's occurrences at that position, IGNORING gaps
     ('-') and ambiguous residues ('X'). Proportions are computed only
     over real (standard) amino acid calls at that position.
  2. Sort the observed amino acids at that position by proportion,
     descending.
  3. The single highest-proportion amino acid is treated as the
     "original" / wild-type residue at that position.
  4. The remaining proportion mass (1 - top_proportion) is halved to get
     a threshold. Walk down the sorted list (excluding the top one) and
     keep adding amino acids to that position's "extra" list,
     accumulating their proportions, for as long as the running total
     does not exceed the threshold. The first amino acid that would push
     the cumulative sum over the threshold is excluded, and no further
     (smaller-proportion) amino acids are considered for that position
     (strict sequential walk, not a best-fit search).

OUTPUT:
  A single-column CSV ("sequence"):
    - the first row is the full consensus sequence, built from the
      top-scoring amino acid at every position
    - one additional row for EVERY (position, extra amino acid) pair
      that passed the 50%-of-remaining threshold. Each such row is the
      consensus sequence with ONLY that one position swapped to the
      alternate amino acid -- every other position stays at its
      top-scoring/original residue.

  Example, for a 12-residue motif where position 1 has two extra
  amino acids (T, M) and position 2 has one extra amino acid (S):

      sequence
      AVASIOCNAIJF
      TVASIOCNAIJF
      MVASIOCNAIJF
      ASASIOCNAIJF

  (A per-position breakdown with labels like "pos 1a"/"pos 1b" is
  still available via the optional --report TSV, for inspection.)

  Optionally, a separate TSV report can also be written, giving the
  full proportion breakdown computed at every position (useful for
  double-checking the thresholding).

NOTES / EDGE CASES:
  - All input sequences must be the same length (i.e. already aligned).
    The script raises an error if they are not.
  - If a position has zero real (non-gap/non-X) residues across every
    input sequence, there's nothing to build a distribution from. That
    position is marked 'X' in the consensus and contributes no variant
    rows. This is logged clearly in the report.
  - Floating-point comparisons use a small epsilon tolerance to avoid
    borderline proportions being excluded/included due to rounding.

INPUT:  a CSV file with a column of sequences (default column name:
        "sequence").

USAGE:
    python filter_sequences.py input.csv output.csv \
        [--seq-col sequence] [--report report.tsv]
"""

import argparse
import csv
import string
import sys
from collections import Counter

GAP_CHARS = {"-", "."}
AMBIGUOUS_CHARS = {"X", "x"}
EPSILON = 1e-9


def read_sequences(csv_path, seq_col):
    """Read sequences from a CSV file."""
    seqs = []
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if seq_col not in reader.fieldnames:
            raise ValueError(
                f"Column '{seq_col}' not found in {csv_path}. "
                f"Available columns: {reader.fieldnames}"
            )
        for row in reader:
            seq = row[seq_col].strip().upper()
            if seq:
                seqs.append(seq)
    if not seqs:
        raise ValueError("No sequences found in input CSV.")
    return seqs


def validate_equal_length(seqs):
    lengths = {len(s) for s in seqs}
    if len(lengths) > 1:
        raise ValueError(
            f"Input sequences are not all the same length (found lengths: "
            f"{sorted(lengths)}). They must already be aligned."
        )
    return lengths.pop()


def summarize_positions(seqs, seq_len):
    """
    Returns:
        position_info: list of dicts, one per position:
                        {'top': aa or None, 'extra': [aa, aa, ...]}
                        'extra' is ordered by descending proportion.
        report_rows:   list of dicts for the optional TSV report.
    """
    position_info = []
    report_rows = []

    for pos in range(seq_len):
        residues = [s[pos] for s in seqs]
        counts = Counter(
            r for r in residues if r not in GAP_CHARS and r not in AMBIGUOUS_CHARS
        )
        total = sum(counts.values())

        if total == 0:
            position_info.append({"top": None, "extra": []})
            report_rows.append(
                {
                    "position": pos + 1,
                    "top": "",
                    "extra": "",
                    "breakdown": "ALL GAP/AMBIGUOUS COLUMN",
                }
            )
            continue

        proportions = sorted(
            ((aa, c / total) for aa, c in counts.items()),
            key=lambda x: -x[1],
        )

        top_aa, top_p = proportions[0]
        remaining = 1.0 - top_p
        threshold = remaining * 0.5

        extra = []
        cumulative = 0.0
        for aa, p in proportions[1:]:
            if cumulative + p <= threshold + EPSILON:
                extra.append(aa)
                cumulative += p
            else:
                break  # sequential walk: stop at first that overshoots

        position_info.append({"top": top_aa, "extra": extra})
        breakdown_str = ", ".join(f"{aa}:{p:.1%}" for aa, p in proportions)
        report_rows.append(
            {
                "position": pos + 1,
                "top": top_aa,
                "extra": "".join(extra),
                "breakdown": breakdown_str,
            }
        )

    return position_info, report_rows


def build_variant_rows(position_info):
    """
    Build the (label, sequence) rows: the consensus/original sequence,
    plus one single-position-swapped variant per extra amino acid.
    """
    consensus = "".join(info["top"] if info["top"] else "X" for info in position_info)

    rows = [("original", consensus)]
    letters = string.ascii_lowercase

    for pos, info in enumerate(position_info):
        for idx, aa in enumerate(info["extra"]):
            suffix = letters[idx] if idx < len(letters) else f"_{idx}"
            label = f"pos {pos + 1}{suffix}"
            variant_seq = consensus[:pos] + aa + consensus[pos + 1:]
            rows.append((label, variant_seq))

    return rows, consensus


def count_matches(seqs, position_info):
    """How many of the original input sequences are fully consistent
    with the derived allowed sets (top + extra) at every position?
    Informational only -- not written to the output CSV."""
    allowed_sets = [
        (None if info["top"] is None else {info["top"], *info["extra"]})
        for info in position_info
    ]
    count = 0
    for seq in seqs:
        ok = True
        for pos, aa in enumerate(seq):
            allowed = allowed_sets[pos]
            if allowed is not None and aa not in allowed:
                ok = False
                break
        if ok:
            count += 1
    return count


def write_variant_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"])
        for _label, seq in rows:
            writer.writerow([seq])


def write_report(path, report_rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["position", "top_aa", "extra_aa", "proportion_breakdown"])
        for row in report_rows:
            writer.writerow([row["position"], row["top"], row["extra"], row["breakdown"]])


def main():
    parser = argparse.ArgumentParser(
        description="Build a degenerate-motif variant table from aligned protein sequences."
    )
    parser.add_argument("input_csv", help="Path to input CSV file")
    parser.add_argument("output_csv", help="Path to write the variant table (CSV)")
    parser.add_argument(
        "--seq-col", default="sequence", help="Name of the sequence column (default: 'sequence')"
    )
    parser.add_argument(
        "--report", default=None, help="Optional path to write a per-position TSV report"
    )
    args = parser.parse_args()

    seqs = read_sequences(args.input_csv, args.seq_col)
    seq_len = validate_equal_length(seqs)

    position_info, report_rows = summarize_positions(seqs, seq_len)
    rows, consensus = build_variant_rows(position_info)

    write_variant_csv(args.output_csv, rows)
    if args.report:
        write_report(args.report, report_rows)

    matches = count_matches(seqs, position_info)

    print(f"Input sequences:        {len(seqs)}", file=sys.stderr)
    print(f"Sequence length:        {seq_len}", file=sys.stderr)
    print(f"Consensus sequence:     {consensus}", file=sys.stderr)
    print(f"Variant rows written:   {len(rows)} (1 original + {len(rows) - 1} variants)", file=sys.stderr)
    print(f"Input seqs matching the derived motif: {matches}/{len(seqs)}", file=sys.stderr)
    print(f"Output written to:      {args.output_csv}", file=sys.stderr)
    if args.report:
        print(f"Report written to:      {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()