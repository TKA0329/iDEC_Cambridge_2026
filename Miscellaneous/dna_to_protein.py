#!/usr/bin/env python3
"""
dna_to_protein.py

Translates DNA sequence(s) to amino acid sequence using the standard
genetic code. Works as:
  1) A one-off translator for a single sequence typed on the command line
  2) A batch CSV translator (same "sequence" column convention as the
     other CAHS pipeline scripts)

USAGE:
    # Single sequence, printed to terminal
    python dna_to_protein.py ATGCGTAAACAACAAGAAGTTGAAGCTGAT

    # Batch mode: CSV in -> CSV out, translates every row's "dna_sequence"
    # column and appends a "protein_sequence" column
    python dna_to_protein.py input.csv output.csv

FRAME / STRAND:
    Translation starts at position 0 of whatever sequence you give it
    (reading frame +1, forward strand). If your insert isn't already
    in-frame from position 0, trim/pad it first, or use --frame /
    --revcomp below.

STOP CODONS:
    By default translation stops at the first in-frame stop codon and
    the stop is NOT included in the output (matches how you'd want a
    CDS translated for expression). Use --readthrough to instead
    translate the whole sequence and mark stops as "*".
"""

import sys
import csv
import argparse

# ---------------------------------------------------------------------------
# Standard genetic code (NCBI translation table 1). Manual table used so
# this script has no dependencies -- gives identical results to
# Bio.Seq.translate(table=1) for standard codons.
# ---------------------------------------------------------------------------
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def clean_dna(seq: str) -> str:
    """Uppercase, strip whitespace, convert U->T (in case of RNA input)."""
    seq = "".join(seq.upper().split())
    return seq.replace("U", "T")


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq))


def translate(seq: str, frame: int = 0, readthrough: bool = False):
    """
    Translates a DNA string starting at the given 0-indexed frame offset
    (0, 1, or 2). Returns (protein_string, warnings_list).
    """
    warnings = []
    seq = seq[frame:]

    trailing = len(seq) % 3
    if trailing != 0:
        warnings.append(
            f"sequence length not a multiple of 3 after frame offset "
            f"({trailing} leftover base{'s' if trailing != 1 else ''} trimmed from the end)"
        )
        seq = seq[: len(seq) - trailing]

    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i : i + 3]
        aa = CODON_TABLE.get(codon)
        if aa is None:
            warnings.append(f"unrecognized codon '{codon}' at position {i} -> 'X'")
            aa = "X"
        if aa == "*" and not readthrough:
            break
        protein.append(aa)

    return "".join(protein), warnings


def run_single(args):
    seq = clean_dna(args.input)
    frame = args.frame
    if args.revcomp:
        seq = reverse_complement(seq)

    protein, warnings = translate(seq, frame=frame, readthrough=args.readthrough)

    print(f"DNA ({len(seq)} nt){' [revcomp]' if args.revcomp else ''}: {seq}")
    print(f"Protein ({len(protein)} aa): {protein}")
    for w in warnings:
        print(f"  WARNING: {w}")


def run_batch(in_path, out_path, frame: int, revcomp: bool, readthrough: bool):
    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames_lower = {name.lower(): name for name in reader.fieldnames}

        dna_col = fieldnames_lower.get("dna_sequence") or fieldnames_lower.get("sequence")
        if dna_col is None:
            print("ERROR: input CSV must have a 'dna_sequence' (or 'sequence') column.")
            sys.exit(1)

        rows_in = list(reader)
        fieldnames = [name for name in reader.fieldnames if name != dna_col]
        fieldnames.append("protein_sequence")
        fieldnames.append("translation_warnings")

    results = []
    n_warned = 0
    for row in rows_in:
        raw = row.pop(dna_col, "")
        seq = clean_dna(raw)
        if revcomp:
            seq = reverse_complement(seq)
        protein, warnings = translate(seq, frame=frame, readthrough=readthrough)
        row["protein_sequence"] = protein
        row["translation_warnings"] = "; ".join(warnings)
        if warnings:
            n_warned += 1
        results.append(row)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Translated {len(results)} sequences -> {out_path}")
    if n_warned:
        print(f"  {n_warned} row(s) had warnings (non-multiple-of-3 length and/or unrecognized codons) -- see 'translation_warnings' column.")


def main():
    parser = argparse.ArgumentParser(
        description="Translate DNA to protein: single sequence or batch CSV mode.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        help="A raw DNA sequence (single mode) OR path to an input CSV (batch mode).",
    )
    parser.add_argument(
        "output_csv",
        nargs="?",
        default=None,
        help="Output CSV path. If given, input is treated as a CSV file (batch mode).",
    )
    parser.add_argument(
        "--frame", type=int, choices=[0, 1, 2], default=0,
        help="0-indexed reading frame offset (0, 1, or 2). Default: 0.",
    )
    parser.add_argument(
        "--revcomp", action="store_true",
        help="Translate the reverse complement instead of the given strand.",
    )
    parser.add_argument(
        "--readthrough", action="store_true",
        help="Translate through stop codons (marks them as '*') instead of stopping at the first one.",
    )
    args = parser.parse_args()

    if args.output_csv:
        run_batch(args.input, args.output_csv, args.frame, args.revcomp, args.readthrough)
    else:
        run_single(args)


if __name__ == "__main__":
    main()