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
    --unique        Sample distinct DNA sequences, no DNA-level repeats (default: on)
    --no-unique     Allow the same DNA sequence to be drawn more than once
                     (--unique fails/warns if the degenerate space is smaller than --num)

PROTEIN-LEVEL DEDUPLICATION (always on, not a flag):
    Two different DNA sequences can translate to the same protein (synonymous
    codons), so DNA-level uniqueness alone does not guarantee unique protein
    output. This script always keeps sampling additional DNA sequences until
    it has --num sequences that are unique AT THE PROTEIN LEVEL (or until the
    reachable protein diversity is exhausted, in which case it warns and
    writes however many unique proteins it actually found). There is no
    --no-protein-dedup escape hatch: the output CSV's "sequence" column is
    always duplicate-free.
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
    seen_proteins = set()

    # Protein-level dedup is ALWAYS on (not a flag): --unique/--no-unique below
    # only control whether the underlying DNA draws are distinct from one
    # another. Even with --unique, two different DNA indices can translate to
    # the same protein (synonymous codons), so we keep drawing DNA until we
    # have --num sequences whose PROTEIN is unique, discarding any draw whose
    # protein has already been seen.
    #
    # Stall breaker: a heavily degenerate motif (e.g. lots of NNK codons) can
    # have a DNA space many orders of magnitude larger than its reachable
    # protein space. Without a circuit-breaker, "keep sampling until unique"
    # would spin through an astronomically large but repetitive DNA space
    # forever, chasing a protein target it can never reach. If this many
    # consecutive draws in a row fail to produce a new unique protein, we
    # conclude the reachable protein diversity is exhausted and stop early.
    STALL_LIMIT = max(20_000, args.num * 20)
    consecutive_no_new = 0
    attempts = 0
    report_interval = max(1, args.num // 20)
    start_time = time.time()

    def report_progress():
        elapsed = time.time() - start_time
        print(
            f"[progress] {len(rows):,}/{args.num:,} unique proteins found "
            f"({len(rows) / args.num:.1%}) | {attempts:,} DNA draws attempted | {elapsed:.1f}s",
            file=sys.stderr,
        )
        sys.stderr.flush()

    if args.unique:
        if args.num > total_space:
            print(
                f"Warning: requested {args.num:,} unique DNA draws but the degenerate "
                f"sequence only encodes {total_space:,}. At most {total_space:,} distinct "
                f"DNA sequences (and therefore at most that many unique proteins) can be produced.",
                file=sys.stderr,
            )
        dna_target = min(args.num, total_space)
        print(f"[run] Drawing distinct DNA sequences (up to {dna_target:,}) from a space of "
              f"{total_space:,}, keeping only new unique proteins...", file=sys.stderr)
        sys.stderr.flush()

        # NOTE: we deliberately do NOT use random.sample(range(total_space), k)
        # here. random.sample() calls len() on its population, and range.__len__
        # raises OverflowError once total_space exceeds a C ssize_t (~9.2e18) --
        # which real degenerate sequences of even moderate length blow past
        # instantly (e.g. 156 nt of mixed N/K codons here works out to ~3.3e25).
        # random.randrange() has no such limit, so we track drawn indices by hand.
        selected_indices = set()
        while len(rows) < args.num and len(selected_indices) < total_space:
            idx = random.randrange(total_space)
            if idx in selected_indices:
                continue
            selected_indices.add(idx)
            attempts += 1
            dna = decode_index(idx, degenerate_seq)
            protein = translate(dna)
            if protein in seen_proteins:
                consecutive_no_new += 1
                if consecutive_no_new >= STALL_LIMIT:
                    print(
                        f"Warning: {STALL_LIMIT:,} consecutive DNA draws in a row produced no new "
                        f"unique protein. Reachable protein diversity appears exhausted; stopping early "
                        f"with {len(rows):,}/{args.num:,} unique proteins found.",
                        file=sys.stderr,
                    )
                    break
                continue
            consecutive_no_new = 0
            seen_proteins.add(protein)
            rows.append((len(rows) + 1, dna, protein))
            if len(rows) % report_interval == 0 or len(rows) == args.num:
                report_progress()

        if len(rows) < args.num and len(selected_indices) >= total_space:
            print(
                f"Warning: the entire degenerate DNA space ({total_space:,} sequences) was drawn "
                f"and only {len(rows):,} unique proteins were found (requested {args.num:,}).",
                file=sys.stderr,
            )
    else:
        print(f"[run] Sampling DNA (with replacement, DNA-level repeats allowed) until "
              f"{args.num:,} unique proteins are found...", file=sys.stderr)
        sys.stderr.flush()
        while len(rows) < args.num:
            dna = sample_dna(degenerate_seq)
            attempts += 1
            protein = translate(dna)
            if protein in seen_proteins:
                consecutive_no_new += 1
                if consecutive_no_new >= STALL_LIMIT:
                    print(
                        f"Warning: {STALL_LIMIT:,} consecutive draws in a row produced no new "
                        f"unique protein. Reachable protein diversity appears exhausted; stopping early "
                        f"with {len(rows):,}/{args.num:,} unique proteins found.",
                        file=sys.stderr,
                    )
                    break
                continue
            consecutive_no_new = 0
            seen_proteins.add(protein)
            rows.append((len(rows) + 1, dna, protein))
            if len(rows) % report_interval == 0 or len(rows) == args.num:
                report_progress()

    if len(rows) < args.num:
        print(
            f"Warning: could only generate {len(rows):,} unique-protein sequences "
            f"(requested {args.num:,}) after {attempts:,} DNA draws.",
            file=sys.stderr,
        )

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"])
        writer.writerows([(protein,) for _seq_id, _dna, protein in rows])

    print(f"Wrote {len(rows)} sequences to {args.output}")


if __name__ == "__main__":
    main()