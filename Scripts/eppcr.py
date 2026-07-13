#!/usr/bin/env python3
"""
Simulate an error-prone PCR mutant library.

Error-prone PCR typically introduces largely random, independent point
substitutions at each nucleotide position. This simulation treats the provided
error rate as the final per-base mutation probability in the amplified product.

For example, 0.6% means 0.006 mutations per nucleotide, or 6 mutations per kb.

The final mutation count for a molecule is drawn from a binomial distribution
over the sequence length, treating every position as an independent site that
may acquire a single substitution.

Optionally, molecule-to-molecule rate variability can be introduced with a
gamma distribution.

Requirements:
    Python 3.9+
    numpy

Install dependency:
    pip install numpy

Example:
    python3 eppcr.py \
        --sequence ATGGCTGCTGCTTAA \
        --count 100000 \
        --error-rate 0.006 \
        --cycles 25 \
        --output mutants.csv \
        --summary mutations.tsv \
        --seed 42
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


VALID_BASES = frozenset("ACGT")
TRANSITION = {
    "A": "G",
    "G": "A",
    "C": "T",
    "T": "C",
}
TRANSVERSIONS = {
    "A": ("C", "T"),
    "G": ("C", "T"),
    "C": ("A", "G"),
    "T": ("A", "G"),
}


@dataclass(frozen=True)
class Mutation:
    position: int       # 1-based position
    reference: str
    alternate: str

    def label(self) -> str:
        return f"{self.reference}{self.position}{self.alternate}"


def clean_sequence(raw_sequence: str) -> str:
    """
    Remove FASTA headers and whitespace, then validate the DNA sequence.
    """
    sequence_lines = [
        line.strip()
        for line in raw_sequence.splitlines()
        if line.strip() and not line.startswith(">")
    ]
    sequence = re.sub(r"\s+", "", "".join(sequence_lines)).upper()

    if not sequence:
        raise ValueError("The input DNA sequence is empty.")

    invalid = sorted(set(sequence) - VALID_BASES)
    if invalid:
        raise ValueError(
            "Only unambiguous A, C, G, and T bases are supported. "
            f"Invalid characters: {', '.join(invalid)}"
        )

    return sequence


def load_sequence(sequence_argument: str | None, input_file: Path | None) -> str:
    if sequence_argument is not None:
        return clean_sequence(sequence_argument)

    if input_file is not None:
        try:
            return clean_sequence(input_file.read_text(encoding="utf-8"))
        except OSError as exc:
            raise ValueError(f"Could not read input file: {exc}") from exc

    raise ValueError("Provide either --sequence or --input.")


def choose_alternate_base(
    reference: str,
    transition_probability: float,
    rng: np.random.Generator,
) -> str:
    """
    Choose a substitution with a specified probability of being a transition.

    Transitions:
        A <-> G
        C <-> T

    Transversions:
        purine <-> pyrimidine
    """
    if rng.random() < transition_probability:
        return TRANSITION[reference]

    alternatives = TRANSVERSIONS[reference]
    return alternatives[int(rng.integers(0, len(alternatives)))]


def normalize_error_rate(error_rate: float) -> float:
    """
    Normalize the user-facing error rate.

    The model uses substitution probabilities in decimal form. To keep the CLI
    ergonomic for common EP-PCR ranges, values in the typical manufacturer
    range of 0.1-100 are interpreted as percentages (for example 0.6 -> 0.006,
    2.0 -> 0.02), while smaller values are treated as decimal probabilities
    (for example 0.006 -> 0.006, 0.02 -> 0.02).
    """
    if error_rate < 0:
        raise ValueError("--error-rate must be non-negative.")

    if 0 <= error_rate <= 0.1:
        return float(error_rate)

    if 0.1 < error_rate <= 100:
        return float(error_rate / 100.0)

    raise ValueError(
        "--error-rate must be between 0 and 100, where values at or below "
        "0.1 are treated as decimal probabilities and larger values are "
        "treated as percentages."
    )


def draw_effective_rate(
    mean_rate: float,
    rate_shape: float | None,
    rng: np.random.Generator,
) -> float:
    """
    Draw a molecule-specific error rate.

    When rate_shape is omitted, every molecule receives the same mean rate.

    When rate_shape is supplied, rates follow a gamma distribution:
        shape = rate_shape
        scale = mean_rate / rate_shape

    Larger shape values produce less molecule-to-molecule variability.
    """
    if rate_shape is None:
        return mean_rate

    return float(rng.gamma(shape=rate_shape, scale=mean_rate / rate_shape))


def mutate_sequence(
    sequence: str,
    per_base_per_cycle_rate: float,
    cycles: int,
    transition_probability: float,
    rate_shape: float | None,
    rng: np.random.Generator,
) -> tuple[str, list[Mutation], float]:
    """
    Generate one mutant sequence.

    Mutation positions are sampled without replacement, so a nucleotide cannot
    be mutated twice in one final output molecule.
    """
    sequence_length = len(sequence)
    effective_rate = draw_effective_rate(
        mean_rate=per_base_per_cycle_rate,
        rate_shape=rate_shape,
        rng=rng,
    )

    mutation_probability = effective_rate
    mutation_count = int(rng.binomial(sequence_length, mutation_probability))

    # The no-repeat-position model cannot exceed the sequence length.
    mutation_count = min(mutation_count, sequence_length)

    if mutation_count == 0:
        return sequence, [], effective_rate

    positions = np.sort(
        rng.choice(sequence_length, size=mutation_count, replace=False)
    )

    mutated = list(sequence)
    mutations: list[Mutation] = []

    for zero_based_position in positions:
        position = int(zero_based_position)
        reference = sequence[position]
        alternate = choose_alternate_base(
            reference=reference,
            transition_probability=transition_probability,
            rng=rng,
        )

        mutated[position] = alternate
        mutations.append(
            Mutation(
                position=position + 1,
                reference=reference,
                alternate=alternate,
            )
        )

    return "".join(mutated), mutations, effective_rate


def simulate_library(
    sequence: str,
    count: int,
    per_base_per_cycle_rate: float,
    cycles: int,
    transition_probability: float,
    rate_shape: float | None,
    allow_wild_type: bool,
    unique_only: bool,
    seed: int | None,
    max_attempts: int,
) -> list[tuple[str, list[Mutation], float]]:
    rng = np.random.default_rng(seed)

    results: list[tuple[str, list[Mutation], float]] = []
    seen_sequences: set[str] = set()
    attempts = 0

    while len(results) < count:
        attempts += 1

        if attempts > max_attempts:
            raise RuntimeError(
                f"Unable to generate {count} qualifying sequences after "
                f"{max_attempts:,} attempts. This usually happens when "
                "--unique-only or --exclude-wild-type is used with a very low "
                "mutation rate. Increase the error rate, cycles, or max attempts."
            )

        mutant, mutations, effective_rate = mutate_sequence(
            sequence=sequence,
            per_base_per_cycle_rate=per_base_per_cycle_rate,
            cycles=cycles,
            transition_probability=transition_probability,
            rate_shape=rate_shape,
            rng=rng,
        )

        if not allow_wild_type and not mutations:
            continue

        if unique_only and mutant in seen_sequences:
            continue

        seen_sequences.add(mutant)
        results.append((mutant, mutations, effective_rate))

    return results


def write_csv(
    output_path: Path,
    results: list[tuple[str, list[Mutation], float]],
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sequence_id",
                "mutation_count",
                "mutations",
                "effective_error_rate_per_base_per_cycle",
                "sequence",
            ]
        )

        for index, (sequence, mutations, effective_rate) in enumerate(
            results, start=1
        ):
            writer.writerow(
                [
                    f"mutant_{index:06d}",
                    len(mutations),
                    ",".join(mutation.label() for mutation in mutations)
                    or "WT",
                    f"{effective_rate:.10g}",
                    sequence,
                ]
            )


def write_summary(
    summary_path: Path,
    results: list[tuple[str, list[Mutation], float]],
) -> None:
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sequence_id",
                "mutation_count",
                "mutations",
                "effective_error_rate_per_base_per_cycle",
                "sequence",
            ]
        )

        for index, (sequence, mutations, effective_rate) in enumerate(
            results, start=1
        ):
            writer.writerow(
                [
                    f"mutant_{index:06d}",
                    len(mutations),
                    ",".join(mutation.label() for mutation in mutations)
                    or "WT",
                    f"{effective_rate:.10g}",
                    sequence,
                ]
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Simulate an error-prone PCR mutant DNA library."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--sequence",
        help="DNA sequence supplied directly on the command line.",
    )
    input_group.add_argument(
        "--input",
        type=Path,
        help="Text or FASTA file containing the DNA sequence.",
    )

    parser.add_argument(
        "--count",
        type=int,
        required=True,
        help="Number of output sequences to generate.",
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        required=True,
        help=(
            "Mean substitution probability per base per PCR cycle. Typical "
            "error-prone PCR values are 0.006-0.02 (0.6-2.0%%). Values "
            "at or below 0.1 are treated as decimal probabilities, while "
            "values above 0.1 are treated as percentages."
        ),
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=1,
        help=(
            "Number of simulated PCR cycles. The mutation rate is interpreted "
            "as the final per-base probability in the amplified product, so "
            "cycle count does not change the expected number of substitutions. "
            "Default: 1."
        ),
    )
    parser.add_argument(
        "--transition-probability",
        type=float,
        default=1.0 / 3.0,
        help=(
            "Probability that a substitution is a transition rather than a "
            "transversion. Default: 0.3333."
        ),
    )
    parser.add_argument(
        "--rate-shape",
        type=float,
        default=None,
        help=(
            "Optional gamma-distribution shape for molecule-to-molecule error "
            "rate variation. Smaller values give more variation. Values such "
            "as 2 to 10 are reasonable exploratory settings."
        ),
    )
    parser.add_argument(
        "--exclude-wild-type",
        action="store_true",
        help="Reject sequences containing zero mutations.",
    )
    parser.add_argument(
        "--unique-only",
        action="store_true",
        help="Require every output DNA sequence to be unique.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible output.",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=1_000_000,
        help="Maximum generation attempts. Default: 1,000,000.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("eppcr_mutants.csv"),
        help="Output CSV filename. Default: eppcr_mutants.csv.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("eppcr_mutations.tsv"),
        help="Output TSV summary filename. Default: eppcr_mutations.tsv. This file is a compact tab-separated summary of the same results.",
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.count <= 0:
        raise ValueError("--count must be greater than zero.")

    args.error_rate = normalize_error_rate(args.error_rate)

    if args.cycles <= 0:
        raise ValueError("--cycles must be greater than zero.")

    if not 0 <= args.transition_probability <= 1:
        raise ValueError("--transition-probability must be between 0 and 1.")

    if args.rate_shape is not None and args.rate_shape <= 0:
        raise ValueError("--rate-shape must be greater than zero.")

    if args.max_attempts <= 0:
        raise ValueError("--max-attempts must be greater than zero.")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        validate_arguments(args)
        sequence = load_sequence(args.sequence, args.input)

        results = simulate_library(
            sequence=sequence,
            count=args.count,
            per_base_per_cycle_rate=args.error_rate,
            cycles=args.cycles,
            transition_probability=args.transition_probability,
            rate_shape=args.rate_shape,
            allow_wild_type=not args.exclude_wild_type,
            unique_only=args.unique_only,
            seed=args.seed,
            max_attempts=args.max_attempts,
        )

        write_csv(args.output, results)
        write_summary(args.summary, results)

    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))

    mutation_counts = np.array(
        [len(mutations) for _, mutations, _ in results],
        dtype=float,
    )

    print(f"Reference length: {len(sequence):,} bp")
    print(f"Sequences generated: {len(results):,}")
    print(f"Mean mutations per sequence: {mutation_counts.mean():.3f}")
    print(f"Median mutations per sequence: {np.median(mutation_counts):.3f}")
    print(f"Minimum mutations: {int(mutation_counts.min())}")
    print(f"Maximum mutations: {int(mutation_counts.max())}")
    print(f"CSV output: {args.output}")
    print(f"Summary output: {args.summary}")


if __name__ == "__main__":
    main()