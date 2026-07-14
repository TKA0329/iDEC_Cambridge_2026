#!/usr/bin/env python3
"""
iupac_sampler.py

Takes a degenerate IUPAC nucleotide sequence and generates random
concrete sample sequences (resolving each ambiguous base to one of
its possible nucleotides), writing them out as a FASTA file.

Usage (interactive):
    python3 iupac_sampler.py

Usage (non-interactive, via arguments):
    python3 iupac_sampler.py --seq RYKSWMBDHVN --num 10 --out samples.fasta
"""

import argparse
import random
import sys

# IUPAC nucleotide ambiguity codes -> possible bases
IUPAC_CODES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def validate_sequence(seq: str) -> str:
    """Uppercase and validate that every character is a recognised IUPAC code."""
    seq = seq.strip().upper()
    invalid = set(seq) - set(IUPAC_CODES.keys())
    if invalid:
        raise ValueError(
            f"Sequence contains invalid IUPAC character(s): {', '.join(sorted(invalid))}"
        )
    if not seq:
        raise ValueError("Sequence cannot be empty.")
    return seq


def generate_sample(seq: str) -> str:
    """Resolve each ambiguous base to a single randomly chosen concrete base."""
    return "".join(random.choice(IUPAC_CODES[base]) for base in seq)


def write_fasta(samples, out_path: str, seq_id_prefix: str = "sample"):
    with open(out_path, "w") as f:
        for i, s in enumerate(samples, start=1):
            f.write(f">{seq_id_prefix}_{i}\n{s}\n")


def get_args():
    parser = argparse.ArgumentParser(
        description="Generate random sample sequences from a degenerate IUPAC sequence."
    )
    parser.add_argument("--seq", help="Degenerate IUPAC sequence (e.g. RYKSWMBDHVN)")
    parser.add_argument("--num", type=int, help="Number of sample sequences to generate")
    parser.add_argument("--out", default="samples.fasta", help="Output FASTA file path")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (optional, for reproducibility)")
    return parser.parse_args()


def main():
    args = get_args()

    if args.seed is not None:
        random.seed(args.seed)

    # Get sequence
    seq = args.seq
    if not seq:
        seq = input("Enter the degenerate IUPAC sequence: ")
    try:
        seq = validate_sequence(seq)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Get number of samples
    num = args.num
    if num is None:
        while True:
            raw = input("How many sample sequences would you like to generate? ")
            try:
                num = int(raw)
                if num <= 0:
                    print("Please enter a positive integer.")
                    continue
                break
            except ValueError:
                print("Please enter a valid integer.")

    # Output path
    out_path = args.out
    if not args.seq:  # only prompt for output name if we were fully interactive
        custom_out = input(f"Output FASTA filename [default: {out_path}]: ").strip()
        if custom_out:
            out_path = custom_out

    samples = [generate_sample(seq) for _ in range(num)]
    write_fasta(samples, out_path)

    print(f"\nDone! {num} sample sequence(s) written to '{out_path}'.")
    print(f"Template sequence: {seq}  (length {len(seq)})")


if __name__ == "__main__":
    main()
