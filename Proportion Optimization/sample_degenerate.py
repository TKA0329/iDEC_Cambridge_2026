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
import time

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


def compute_space(degenerate_seq: str) -> int:
    """Total number of concrete DNA sequences the degenerate sequence encodes
    (product of the number of IUPAC options at each position)."""
    total = 1
    for base in degenerate_seq:
        total *= len(IUPAC_CODES[base])
    return total


def decode_index(index: int, degenerate_seq: str) -> str:
    """Deterministically map an integer in [0, compute_space(seq)) to one
    concrete DNA sequence (mixed-radix decomposition, one digit per
    position). Every index maps to exactly one sequence and vice versa, so
    sampling distinct indices guarantees distinct sequences with no
    duplicate-checking needed."""
    bases = []
    for base in reversed(degenerate_seq):
        options = IUPAC_CODES[base]
        n_opts = len(options)
        index, digit = divmod(index, n_opts)
        bases.append(options[digit])
    return "".join(reversed(bases))


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
            raw_lines = f.readlines()
        # Tolerate FASTA-style input (e.g. the --out-seq file written by
        # reverse_translation_V5_StopSafe.py): strip '>' header lines and
        # blank lines, and join any wrapped sequence lines. If more than
        # one FASTA record is present (multi-bin output), only the first
        # is used -- this script only supports a single degenerate
        # sequence at a time.
        seq_lines = []
        record_count = 0
        for line in raw_lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_count += 1
                if record_count > 1:
                    break
                continue
            seq_lines.append(line)
        if record_count > 1:
            print(
                f"Warning: --seq-file contains {record_count} FASTA records "
                f"(e.g. multiple RT bins); only the first record's sequence "
                f"is being used.",
                file=sys.stderr,
            )
        raw_seq = "".join(seq_lines) if seq_lines else "".join(raw_lines)
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

    total_space = compute_space(degenerate_seq)
    print(f"[info] Degenerate sequence length: {len(degenerate_seq)}", file=sys.stderr)
    print(f"[info] Total combinatorial space: {total_space:,} possible sequences", file=sys.stderr)
    sys.stderr.flush()

    rows = []

    if args.unique:
        k = min(args.num, total_space)
        if args.num > total_space:
            print(
                f"Warning: requested {args.num:,} unique sequences but the degenerate "
                f"sequence only encodes {total_space:,}. Generating all {total_space:,}.",
                file=sys.stderr,
            )

        print(f"[run] Sampling {k:,} unique indices from a space of {total_space:,}...", file=sys.stderr)
        sys.stderr.flush()
        t0 = time.time()
        # Sample distinct integers, not distinct strings: this de-duplicates
        # via a cheap integer set, so it stays fast even when k is close to
        # total_space (no coupon-collector blowup on long DNA strings).
        #
        # NOTE: we deliberately do NOT use random.sample(range(total_space), k)
        # here. random.sample() calls len() on its population, and range.__len__
        # raises OverflowError once total_space exceeds a C ssize_t (~9.2e18) --
        # which real degenerate sequences of even moderate length blow past
        # instantly (e.g. 156 nt of mixed N/K codons here works out to ~3.3e25).
        # random.randrange() has no such limit: it handles arbitrarily large
        # Python ints directly via getrandbits(), so we build the distinct set
        # by hand instead.
        selected = set()
        while len(selected) < k:
            selected.add(random.randrange(total_space))
        indices = list(selected)
        print(f"[run] Index sampling done in {time.time() - t0:.2f}s. Decoding + translating...", file=sys.stderr)
        sys.stderr.flush()

        report_interval = max(1, k // 20)
        start_time = time.time()
        for i, idx in enumerate(indices, start=1):
            dna = decode_index(idx, degenerate_seq)
            protein = translate(dna)
            rows.append((i, dna, protein))
            if i % report_interval == 0 or i == k:
                elapsed = time.time() - start_time
                print(f"[progress] {i:,}/{k:,} decoded ({i / k:.1%}) | {elapsed:.1f}s", file=sys.stderr)
                sys.stderr.flush()
    else:
        report_interval = max(1, args.num // 20)
        start_time = time.time()
        for i in range(1, args.num + 1):
            dna = sample_dna(degenerate_seq)
            protein = translate(dna)
            rows.append((i, dna, protein))
            if i % report_interval == 0 or i == args.num:
                elapsed = time.time() - start_time
                print(f"[progress] {i:,}/{args.num:,} sampled ({i / args.num:.1%}) | {elapsed:.1f}s", file=sys.stderr)
                sys.stderr.flush()

    if args.unique and len(rows) < args.num:
        print(
            f"Warning: could only generate {len(rows)} unique sequences "
            f"(requested {args.num}). The degenerate sequence's total possible "
            f"space may be smaller than requested.",
            file=sys.stderr,
        )

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"])
        writer.writerows([(protein,) for _seq_id, _dna, protein in rows])

    print(f"Wrote {len(rows)} sequences to {args.output}")


if __name__ == "__main__":
    main()