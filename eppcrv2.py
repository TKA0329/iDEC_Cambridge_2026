#!/usr/bin/env python3
"""
Simulate an error-prone PCR mutant library.

Error-prone PCR typically introduces largely random, independent point
substitutions at each nucleotide position. This simulation treats the
provided error rate as the PER-CYCLE per-base mutation probability, and
compounds it across --cycles to get the final per-base probability in the
amplified product:

    final_per_base_probability = 1 - (1 - per_cycle_rate) ** cycles

For example, a per-cycle rate of 0.6% (0.006) over 25 cycles gives a final
per-base probability of ~14.1%, i.e. roughly 141 mutations per kb -- NOT
6 mutations per kb. If you already have a "final" rate in hand (e.g. from a
kit datasheet reporting mutations/kb after N cycles) either set --cycles 1
so no compounding occurs, or back-calculate the per-cycle rate first.

The final mutation count for a molecule is drawn from a binomial distribution
over the sequence length, treating every position as an independent site that
may acquire a single substitution (positions are drawn without replacement,
so no site is hit twice in a single output molecule).

Optionally, molecule-to-molecule rate variability can be introduced with a
gamma distribution, and/or a custom substitution spectrum (transition bias
per reference base) and per-position hotspot weighting can be supplied to
better approximate a specific polymerase's known error bias.

Requirements:
    Python 3.9+
    numpy

Install dependency:
    pip install numpy

Example:
    python3 eppcr.py \
        --sequence ATGGCTGCTGCTTAA \
        --count 100000 \
        --error-rate 0.6 --units percent \
        --cycles 25 \
        --output mutants.csv \
        --summary mutations.tsv \
        --seed 42
"""

from __future__ import annotations

import argparse
import csv
import json
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

    Raises if more than one FASTA record (i.e. more than one header line) is
    present, since silently concatenating multiple records into one sequence
    would be a silent correctness bug rather than a helpful convenience.
    """
    header_count = sum(
        1 for line in raw_sequence.splitlines() if line.strip().startswith(">")
    )
    if header_count > 1:
        raise ValueError(
            f"Input contains {header_count} FASTA records, but this tool "
            "simulates mutants of a single reference sequence. Please "
            "provide a file with exactly one FASTA record (or a single "
            "raw sequence)."
        )

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


def load_bias_matrix(bias_matrix_file: Path | None) -> dict[str, dict[str, float]] | None:
    """
    Load an optional per-reference-base substitution spectrum.

    The JSON file must map each of A/C/G/T to a mapping of the other three
    bases to non-negative relative weights, e.g.:

        {
          "A": {"G": 4.0, "C": 1.0, "T": 1.0},
          "G": {"A": 4.0, "C": 1.0, "T": 1.0},
          "C": {"T": 4.0, "A": 1.0, "G": 1.0},
          "T": {"C": 4.0, "A": 1.0, "G": 1.0}
        }

    Weights are normalized to probabilities per reference base. This lets you
    encode a specific polymerase's known transition/transversion bias (e.g.
    Taq's well-documented preference for particular transitions) instead of
    the single global --transition-probability figure.
    """
    if bias_matrix_file is None:
        return None

    try:
        raw = json.loads(bias_matrix_file.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError(f"Could not read --bias-matrix file: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"--bias-matrix file is not valid JSON: {exc}") from exc

    if set(raw.keys()) != VALID_BASES:
        raise ValueError(
            "--bias-matrix JSON must have exactly the top-level keys A, C, G, T."
        )

    normalized: dict[str, dict[str, float]] = {}
    for reference, alt_weights in raw.items():
        expected_alts = VALID_BASES - {reference}
        if set(alt_weights.keys()) != expected_alts:
            raise ValueError(
                f"--bias-matrix entry for '{reference}' must have exactly the "
                f"keys {sorted(expected_alts)}."
            )

        weights = {alt: float(weight) for alt, weight in alt_weights.items()}
        if any(weight < 0 for weight in weights.values()):
            raise ValueError("--bias-matrix weights must be non-negative.")

        total = sum(weights.values())
        if total <= 0:
            raise ValueError(
                f"--bias-matrix weights for '{reference}' must sum to a "
                "positive value."
            )

        normalized[reference] = {alt: weight / total for alt, weight in weights.items()}

    return normalized


def load_hotspot_weights(
    hotspot_file: Path | None,
    sequence_length: int,
) -> np.ndarray | None:
    """
    Load optional per-position relative mutation-rate multipliers.

    File format: one "position<TAB or space>weight" pair per line, 1-based
    positions. Unlisted positions default to weight 1.0 (baseline rate).
    Weights are relative multipliers, not probabilities -- a weight of 3.0
    means that position is roughly 3x as likely to be selected as a
    baseline-weight position, given the same total mutation count.
    """
    if hotspot_file is None:
        return None

    weights = np.ones(sequence_length, dtype=float)

    try:
        lines = hotspot_file.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise ValueError(f"Could not read --hotspot-weights file: {exc}") from exc

    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        parts = stripped.replace(",", " ").split()
        if len(parts) != 2:
            raise ValueError(
                f"--hotspot-weights line {line_number} must have exactly two "
                f"fields (position, weight): {line!r}"
            )

        try:
            position = int(parts[0])
            weight = float(parts[1])
        except ValueError as exc:
            raise ValueError(
                f"--hotspot-weights line {line_number} could not be parsed: {exc}"
            ) from exc

        if not 1 <= position <= sequence_length:
            raise ValueError(
                f"--hotspot-weights line {line_number} has position {position}, "
                f"which is out of range for a sequence of length {sequence_length}."
            )
        if weight < 0:
            raise ValueError(
                f"--hotspot-weights line {line_number} has a negative weight."
            )

        weights[position - 1] = weight

    return weights


def choose_alternate_base(
    reference: str,
    transition_probability: float,
    bias_matrix: dict[str, dict[str, float]] | None,
    rng: np.random.Generator,
) -> str:
    """
    Choose a substitution for the given reference base.

    If a bias_matrix was supplied, sample the alternate base directly from
    its normalized per-reference-base distribution. Otherwise fall back to
    the global transition_probability model:

        Transitions:
            A <-> G
            C <-> T

        Transversions:
            purine <-> pyrimidine
    """
    if bias_matrix is not None:
        alt_weights = bias_matrix[reference]
        alternates = list(alt_weights.keys())
        probabilities = list(alt_weights.values())
        return str(rng.choice(alternates, p=probabilities))

    if rng.random() < transition_probability:
        return TRANSITION[reference]

    alternatives = TRANSVERSIONS[reference]
    return alternatives[int(rng.integers(0, len(alternatives)))]


def normalize_error_rate(error_rate: float, units: str) -> float:
    """
    Normalize the user-facing per-cycle error rate to a decimal probability.

    units == "percent":  error_rate is a percentage, e.g. 0.6 -> 0.006
    units == "fraction":  error_rate is already a decimal probability, e.g.
                          0.006 -> 0.006
    """
    if error_rate < 0:
        raise ValueError("--error-rate must be non-negative.")

    if units == "percent":
        if error_rate > 100:
            raise ValueError("--error-rate as a percentage must be <= 100.")
        return float(error_rate / 100.0)

    # units == "fraction"
    if error_rate > 1:
        raise ValueError("--error-rate as a fraction must be <= 1.")
    return float(error_rate)


def compound_across_cycles(per_cycle_rate: float, cycles: int) -> float:
    """
    Compound a per-cycle per-base substitution probability across cycles,
    assuming independence between cycles:

        final_probability = 1 - (1 - per_cycle_rate) ** cycles
    """
    return float(1.0 - (1.0 - per_cycle_rate) ** cycles)


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

    Larger shape values produce less molecule-to-molecule variability. The
    result is clipped to [0, 1] since it is used as a binomial probability.
    """
    if rate_shape is None:
        return mean_rate

    drawn = float(rng.gamma(shape=rate_shape, scale=mean_rate / rate_shape))
    return min(drawn, 1.0)


def mutate_sequence(
    sequence: str,
    final_per_base_rate: float,
    transition_probability: float,
    bias_matrix: dict[str, dict[str, float]] | None,
    rate_shape: float | None,
    hotspot_weights: np.ndarray | None,
    rng: np.random.Generator,
) -> tuple[str, list[Mutation], float]:
    """
    Generate one mutant sequence.

    Mutation positions are sampled without replacement, so a nucleotide cannot
    be mutated twice in one final output molecule. If hotspot_weights are
    supplied, positions are sampled proportionally to those weights instead
    of uniformly.
    """
    sequence_length = len(sequence)
    effective_rate = draw_effective_rate(
        mean_rate=final_per_base_rate,
        rate_shape=rate_shape,
        rng=rng,
    )

    mutation_count = int(rng.binomial(sequence_length, effective_rate))

    if mutation_count == 0:
        return sequence, [], effective_rate

    if hotspot_weights is None:
        positions = np.sort(
            rng.choice(sequence_length, size=mutation_count, replace=False)
        )
    else:
        selection_probabilities = hotspot_weights / hotspot_weights.sum()
        positions = np.sort(
            rng.choice(
                sequence_length,
                size=mutation_count,
                replace=False,
                p=selection_probabilities,
            )
        )

    mutated = list(sequence)
    mutations: list[Mutation] = []

    for zero_based_position in positions:
        position = int(zero_based_position)
        reference = sequence[position]
        alternate = choose_alternate_base(
            reference=reference,
            transition_probability=transition_probability,
            bias_matrix=bias_matrix,
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
    final_per_base_rate: float,
    transition_probability: float,
    bias_matrix: dict[str, dict[str, float]] | None,
    rate_shape: float | None,
    hotspot_weights: np.ndarray | None,
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
            final_per_base_rate=final_per_base_rate,
            transition_probability=transition_probability,
            bias_matrix=bias_matrix,
            rate_shape=rate_shape,
            hotspot_weights=hotspot_weights,
            rng=rng,
        )

        if not allow_wild_type and not mutations:
            continue

        if unique_only and mutant in seen_sequences:
            continue

        seen_sequences.add(mutant)
        results.append((mutant, mutations, effective_rate))

    return results


def write_table(
    output_path: Path,
    results: list[tuple[str, list[Mutation], float]],
    delimiter: str,
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter=delimiter)
        writer.writerow(
            [
                "sequence_id",
                "mutation_count",
                "mutations",
                "effective_error_rate_per_base",
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
        help="Text or single-record FASTA file containing the DNA sequence.",
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
            "PER-CYCLE mean substitution probability per base. Typical "
            "error-prone PCR values are 0.2-2.0 percent per cycle. "
            "Interpreted according to --units."
        ),
    )
    parser.add_argument(
        "--units",
        choices=["percent", "fraction"],
        default="percent",
        help=(
            "How to interpret --error-rate. 'percent' (default): e.g. 0.6 "
            "means 0.6%%. 'fraction': e.g. 0.006 means a decimal probability "
            "of 0.006. There is no auto-detection -- pick the one matching "
            "your input to avoid a 100x unit error."
        ),
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=1,
        help=(
            "Number of simulated PCR cycles. --error-rate is treated as the "
            "PER-CYCLE probability and compounded across cycles as "
            "1 - (1 - rate) ** cycles to get the final per-base probability "
            "in the amplified product. Use --cycles 1 if your --error-rate "
            "already represents a final, post-amplification rate. Default: 1."
        ),
    )
    parser.add_argument(
        "--transition-probability",
        type=float,
        default=1.0 / 3.0,
        help=(
            "Probability that a substitution is a transition rather than a "
            "transversion. Ignored if --bias-matrix is supplied. "
            "Default: 0.3333 (no bias)."
        ),
    )
    parser.add_argument(
        "--bias-matrix",
        type=Path,
        default=None,
        help=(
            "Optional JSON file giving a full per-reference-base substitution "
            "spectrum (relative weights for each possible alternate base), to "
            "model a specific polymerase's known error bias instead of the "
            "single global --transition-probability. Overrides "
            "--transition-probability when supplied. See module docstring "
            "for format."
        ),
    )
    parser.add_argument(
        "--hotspot-weights",
        type=Path,
        default=None,
        help=(
            "Optional file of 'position weight' lines (1-based) giving "
            "relative per-position mutation-rate multipliers, to model known "
            "sequence-context hotspots/coldspots instead of uniform site "
            "selection. Unlisted positions default to weight 1.0."
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
        help=(
            "Output TSV summary filename. Default: eppcr_mutations.tsv. "
            "This file has the same content as --output, tab-separated."
        ),
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    if args.count <= 0:
        raise ValueError("--count must be greater than zero.")

    args.error_rate = normalize_error_rate(args.error_rate, args.units)

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
        bias_matrix = load_bias_matrix(args.bias_matrix)
        hotspot_weights = load_hotspot_weights(args.hotspot_weights, len(sequence))

        final_per_base_rate = compound_across_cycles(args.error_rate, args.cycles)
        if final_per_base_rate > 0.5:
            parser.error(
                f"The compounded per-base mutation probability is "
                f"{final_per_base_rate:.3f} (>50%), which means your sequence "
                "would statistically be scrambled beyond recognition. Check "
                "your --error-rate, --units, and --cycles."
            )

        results = simulate_library(
            sequence=sequence,
            count=args.count,
            final_per_base_rate=final_per_base_rate,
            transition_probability=args.transition_probability,
            bias_matrix=bias_matrix,
            rate_shape=args.rate_shape,
            hotspot_weights=hotspot_weights,
            allow_wild_type=not args.exclude_wild_type,
            unique_only=args.unique_only,
            seed=args.seed,
            max_attempts=args.max_attempts,
        )

        write_table(args.output, results, delimiter=",")
        write_table(args.summary, results, delimiter="\t")

    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))

    mutation_counts = np.array(
        [len(mutations) for _, mutations, _ in results],
        dtype=float,
    )

    print(f"Reference length: {len(sequence):,} bp")
    print(f"Per-cycle per-base rate: {args.error_rate:.6g} ({args.error_rate * 100:.4g}%)")
    print(f"Cycles: {args.cycles}")
    print(f"Final compounded per-base rate: {final_per_base_rate:.6g} ({final_per_base_rate * 100:.4g}%)")
    print(f"Sequences generated: {len(results):,}")
    print(f"Mean mutations per sequence: {mutation_counts.mean():.3f}")
    print(f"Median mutations per sequence: {np.median(mutation_counts):.3f}")
    print(f"Minimum mutations: {int(mutation_counts.min())}")
    print(f"Maximum mutations: {int(mutation_counts.max())}")
    print(f"CSV output: {args.output}")
    print(f"Summary output: {args.summary}")


if __name__ == "__main__":
    main()
