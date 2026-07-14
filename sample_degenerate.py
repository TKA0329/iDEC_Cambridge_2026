#!/usr/bin/env python3
"""
sample_degenerate.py

Takes a degenerate DNA sequence (IUPAC ambiguity codes), randomly samples a
user-specified number of concrete DNA sequences from it, translates each to
protein, and writes the results to a CSV file.

Usage:
    python sample_degenerate.py --sequence "ATGNNKGGCTGA" --num 100 --output samples.csv
    python sample_degenerate.py --seq-file degenerate.txt --num 500 --output samples.csv
    python sample_degenerate.py --sequence "ATGNNKGGC" --num 50 --output out.csv --seed 7
    python sample_degenerate.py --sequence "ATGNNKGGC" --num 50 --output out.csv --no-unique

Options:
    --sequence      Degenerate DNA sequence as a string (IUPAC codes allowed)
    --seq-file      Path to a text file containing the degenerate DNA sequence
                     (use this instead of --sequence)
    --num           Number of sequences to sample (required)
    --output        Output CSV path (default: sampled_sequences.csv)
    --seed          Random seed for reproducibility (default: 42)
    --unique        Ensure all sampled DNA sequences are unique (default: on)
    --no-unique     Allow duplicate sampled DNA sequences
                     (--unique fails/warns if the degenerate space is smaller than --num)
"""

import argparse
import csv
import random
import sys

# IUPAC nucleotide ambiguity codes -> possible bases
IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}

# Standard genetic code (DNA codons -> amino acid, '*' = stop)
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


def validate_sequence(seq: str) -> str:
    seq = seq.strip().upper().replace(" ", "").replace("\n", "")
    invalid = set(seq) - set(IUPAC_CODES.keys())
    if invalid:
        raise ValueError(
            f"Sequence contains invalid characters: {sorted(invalid)}. "
            f"Allowed IUPAC codes: {sorted(IUPAC_CODES.keys())}"
        )
    return seq


def sample_dna(degenerate_seq: str) -> str:
    """Resolve each IUPAC code in the sequence to one randomly chosen base."""
    return "".join(random.choice(IUPAC_CODES[base]) for base in degenerate_seq)


def translate(dna_seq: str) -> str:
    """Translate a DNA sequence to protein. Stops at the first stop codon
    (stop codon itself is not included). Trailing incomplete codon (if the
    sequence length isn't a multiple of 3) is ignored."""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i + 3]
        aa = CODON_TABLE[codon]
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def main():
    parser = argparse.ArgumentParser(
        description="Sample sequences from a degenerate DNA sequence and translate to protein."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sequence", type=str, help="Degenerate DNA sequence string")
    group.add_argument("--seq-file", type=str, help="Path to a file containing the degenerate DNA sequence")

    parser.add_argument("--num", type=int, required=True, help="Number of sequences to sample")
    parser.add_argument("--output", type=str, default="sampled_sequences.csv", help="Output CSV path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)")
    parser.add_argument("--unique", dest="unique", action="store_true", default=True,
                         help="Ensure all sampled DNA sequences are unique (default: on)")
    parser.add_argument("--no-unique", dest="unique", action="store_false",
                         help="Allow duplicate sampled DNA sequences")
    args = parser.parse_args()

    random.seed(args.seed)

    if args.seq_file:
        with open(args.seq_file) as f:
            raw_seq = f.read()
    else:
        raw_seq = args.sequence

    try:
        degenerate_seq = validate_sequence(raw_seq)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if len(degenerate_seq) % 3 != 0:
        print(
            f"Warning: sequence length ({len(degenerate_seq)}) is not a multiple of 3. "
            f"Trailing {len(degenerate_seq) % 3} base(s) will be ignored during translation.",
            file=sys.stderr,
        )

    seen = set()
    rows = []
    max_attempts = args.num * 100 if args.unique else args.num  # safety cap for --unique
    attempts = 0

    while len(rows) < args.num and attempts < max_attempts:
        attempts += 1
        dna = sample_dna(degenerate_seq)
        if args.unique:
            if dna in seen:
                continue
            seen.add(dna)
        protein = translate(dna)
        rows.append((len(rows) + 1, dna, protein))

    if args.unique and len(rows) < args.num:
        print(
            f"Warning: could only generate {len(rows)} unique sequences "
            f"(requested {args.num}). The degenerate sequence's total possible "
            f"space may be smaller than requested.",
            file=sys.stderr,
        )

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequences"])
        writer.writerows([(protein,) for _seq_id, _dna, protein in rows])

    print(f"Wrote {len(rows)} sequences to {args.output}")


if __name__ == "__main__":
    main()