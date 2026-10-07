#!/usr/bin/env python3
"""
motif1_filter_pipeline.py

Streaming scorer + filter for very large mutant libraries (built/tested for
10M-row, 1GB+ CSVs). Scores each sequence exactly like motif1_helix_scorer.py
(Pace & Scholtz helix propensity, Eisenberg hydrophobic moment, helical face
occupancy, salt-bridge count, Pro/Gly hard filters) and writes ONLY the rows
that pass your thresholds -- nothing is held in memory across rows, so this
scales to files far bigger than your RAM.

KEY DIFFERENCE FROM motif1_helix_scorer.py:
- Old script: read all rows into a list, scored them, wrote all rows out.
  Fine for hundreds of rows; at 10M rows this can use many GB of RAM.
- This script: reads one row, scores it, checks filters, writes it (or
  doesn't), and immediately moves on. Memory use stays flat regardless of
  file size.

USAGE:
    python cahs_analyzer_03.py input.csv output.csv

INPUT CSV:
    Must have a column named "sequence" (case-insensitive). An optional
    "id" column is carried through.

OUTPUT CSV:
    Only rows that PASS all active filters, with all scored columns.
    A summary of how many rows were excluded by each criterion is printed
    to the terminal at the end (nothing extra is written to disk, so this
    won't blow up your output size).
"""

import sys
import csv
import math

# ---------------------------------------------------------------------------
# CONSTRUCT NUMBERING (same as before)
# ---------------------------------------------------------------------------
CORE_START_RESNUM = 121
MOTIF1_START_RESNUM = 124
MOTIF1_END_RESNUM = 142

# ---------------------------------------------------------------------------
# FILTER THRESHOLDS
# Set to None to disable any individual bound.
# ---------------------------------------------------------------------------
PS_SUM_MAX = 6.68        # exclude if pace_scholtz_sum  > this (higher = less helix-favorable)
PS_MEAN_MAX = 0.352      # exclude if pace_scholtz_mean > this
PS_MEAN_MIN = None       # OFF by default. Set e.g. 0.10 to enforce an "over-rigidity" floor
                          # (reject mutants that are FAR more helix-favorable than WT).

MOMENT_MIN = 0.263       # exclude if hydrophobic_moment < this
MOMENT_MAX = 0.526        # OFF by default. Set e.g. 0.40 (~1.5x WT) to cap extreme amphipathicity.

FACE_OCCUPANCY_MIN = 0.6  # exclude if helical_face_occupancy < this

EXCLUDE_HELIX_BREAKERS = False  # hard filter: drop any mutant with Pro/Gly in motif 1
                                 # (matches the earlier agreed hard-filter rule; set False to disable)

# ---------------------------------------------------------------------------
# Skip whole-sequence Biopython baseline (GRAVY/instability/etc.) for speed.
# These aren't part of your filter criteria -- set True if you want them
# back in the output (costs extra time across 10M rows).
# ---------------------------------------------------------------------------
COMPUTE_WHOLE_SEQUENCE_BASELINE = False

PROGRESS_EVERY = 200_000  # print a progress line every N input rows

# ---------------------------------------------------------------------------
# Pace & Scholtz (1998) helix propensity scale, kcal/mol.
# ---------------------------------------------------------------------------
PACE_SCHOLTZ = {
    "A": 0.00, "R": 0.21, "N": 0.65, "D": 0.69, "C": 0.68,
    "Q": 0.39, "E": 0.40, "G": 1.00, "H": 0.61, "I": 0.41,
    "L": 0.21, "K": 0.26, "M": 0.24, "F": 0.54, "P": 3.16,
    "S": 0.50, "T": 0.66, "W": 0.49, "Y": 0.53, "V": 0.61,
}

# ---------------------------------------------------------------------------
# Eisenberg consensus hydrophobicity scale.
# ---------------------------------------------------------------------------
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

DELTA_DEG = 100.0

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False


def clean_sequence(seq: str) -> str:
    return "".join(seq.upper().split())


def extract_motif1(full_seq: str) -> str:
    start_idx = MOTIF1_START_RESNUM - CORE_START_RESNUM
    end_idx = MOTIF1_END_RESNUM - CORE_START_RESNUM + 1
    # print(full_seq[start_idx:end_idx])
    return full_seq[start_idx:end_idx]


def pace_scholtz_scores(seq: str):
    vals = [PACE_SCHOLTZ.get(aa) for aa in seq]
    valid = [v for v in vals if v is not None]
    if not valid:
        return None, None
    return sum(valid), sum(valid) / len(valid)


def hydrophobic_moment_and_face(seq: str, delta_deg: float = DELTA_DEG):
    delta_rad = math.radians(delta_deg)
    sin_sum, cos_sum = 0.0, 0.0
    angles, hydrophobic_flags = [], []

    for i, aa in enumerate(seq):
        h = EISENBERG.get(aa)
        if h is None:
            continue
        angle = i * delta_rad
        sin_sum += h * math.sin(angle)
        cos_sum += h * math.cos(angle)
        angles.append(angle)
        hydrophobic_flags.append(h > 0)

    n = len(seq)
    if n == 0:
        return None, None

    moment = math.sqrt(sin_sum ** 2 + cos_sum ** 2) / n
    moment_vector_angle = math.atan2(sin_sum, cos_sum)

    if not any(hydrophobic_flags):
        return moment, 0.0

    on_face = total_hydrophobic = 0
    for angle, is_hydrophobic in zip(angles, hydrophobic_flags):
        if not is_hydrophobic:
            continue
        total_hydrophobic += 1
        diff = abs(math.degrees(angle - moment_vector_angle)) % 360
        diff = min(diff, 360 - diff)
        if diff <= 60:
            on_face += 1

    face_score = on_face / total_hydrophobic if total_hydrophobic else None
    return moment, face_score


def salt_bridge_score(seq: str) -> int:
    acidic, basic = {"E", "D"}, {"K", "R"}
    count = 0
    n = len(seq)
    for i in range(n):
        for offset in (3, 4):
            j = i + offset
            if j >= n:
                continue
            if (seq[i] in acidic and seq[j] in basic) or (seq[i] in basic and seq[j] in acidic):
                count += 1
    return count


def helix_breaker_flags(seq: str):
    n_pro, n_gly = seq.count("P"), seq.count("G")
    return n_pro, n_gly, (n_pro > 0 or n_gly > 0)


def biopython_baseline(seq: str):
    if not HAVE_BIOPYTHON:
        return {}
    try:
        pa = ProteinAnalysis(seq)
        aa_percent = pa.get_amino_acids_percent()
        aliphatic_index = (
            aa_percent.get("A", 0) * 100
            + 2.9 * aa_percent.get("V", 0) * 100
            + 3.9 * (aa_percent.get("I", 0) + aa_percent.get("L", 0)) * 100
        )
        return {
            "gravy": pa.gravy(),
            "instability_index": pa.instability_index(),
            "aliphatic_index": aliphatic_index,
            "net_charge_pH7": pa.charge_at_pH(7.0),
            "isoelectric_point": pa.isoelectric_point(),
        }
    except Exception:
        return {}


def score_sequence(full_seq: str):
    full_seq = clean_sequence(full_seq)
    motif1_seq = extract_motif1(full_seq)

    ps_sum, ps_mean = pace_scholtz_scores(motif1_seq)
    moment, face_score = hydrophobic_moment_and_face(motif1_seq)
    sb_score = salt_bridge_score(motif1_seq)
    n_pro, n_gly, breaker_flag = helix_breaker_flags(motif1_seq)

    row = {
        "motif1_sequence": motif1_seq,
        "full_length": len(full_seq),
        "pace_scholtz_sum": round(ps_sum, 3) if ps_sum is not None else None,
        "pace_scholtz_mean": round(ps_mean, 3) if ps_mean is not None else None,
        "hydrophobic_moment": round(moment, 3) if moment is not None else None,
        "helical_face_occupancy": round(face_score, 3) if face_score is not None else None,
        "salt_bridge_count": sb_score,
        "num_proline": n_pro,
        "num_glycine": n_gly,
        "helix_breaker_flag": breaker_flag,
    }
    if COMPUTE_WHOLE_SEQUENCE_BASELINE:
        row.update(biopython_baseline(full_seq))
    return row


def failing_reason(row):
    """Returns the first filter this row fails, or None if it passes all."""
    if EXCLUDE_HELIX_BREAKERS and row["helix_breaker_flag"]:
        return "helix_breaker"
    if PS_SUM_MAX is not None and row["pace_scholtz_sum"] is not None and row["pace_scholtz_sum"] > PS_SUM_MAX:
        return "pace_scholtz_sum_too_high"
    if PS_MEAN_MAX is not None and row["pace_scholtz_mean"] is not None and row["pace_scholtz_mean"] > PS_MEAN_MAX:
        return "pace_scholtz_mean_too_high"
    if PS_MEAN_MIN is not None and row["pace_scholtz_mean"] is not None and row["pace_scholtz_mean"] < PS_MEAN_MIN:
        return "pace_scholtz_mean_too_low_overrigid"
    if MOMENT_MIN is not None and row["hydrophobic_moment"] is not None and row["hydrophobic_moment"] < MOMENT_MIN:
        return "moment_too_low"
    if MOMENT_MAX is not None and row["hydrophobic_moment"] is not None and row["hydrophobic_moment"] > MOMENT_MAX:
        return "moment_too_high_extreme"
    if FACE_OCCUPANCY_MIN is not None and row["helical_face_occupancy"] is not None and row["helical_face_occupancy"] < FACE_OCCUPANCY_MIN:
        return "face_occupancy_too_low"
    return None


def main():
    if len(sys.argv) != 3:
        print("Usage: python motif1_filter_pipeline.py input.csv output.csv")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    n_total = 0
    n_passed = 0
    n_empty = 0
    exclusion_counts = {}

    out_f = open(out_path, "w", newline="")
    writer = None  # created once we know the fieldnames

    with open(in_path, newline="") as in_f:
        reader = csv.DictReader(in_f)
        fieldnames_lower = {name.lower(): name for name in reader.fieldnames}
        if "sequence" not in fieldnames_lower:
            print("ERROR: input CSV must have a 'sequence' column.")
            sys.exit(1)
        seq_col = fieldnames_lower["sequence"]
        id_col = fieldnames_lower.get("id")

        for i, in_row in enumerate(reader, start=1):
            n_total += 1
            if n_total % PROGRESS_EVERY == 0:
                print(f"  ...processed {n_total:,} rows, {n_passed:,} passed so far")

            raw_seq = (in_row.get(seq_col) or "").strip()
            if not raw_seq:
                n_empty += 1
                continue

            seq_id = in_row.get(id_col, f"seq_{i}") if id_col else f"seq_{i}"
            scored = score_sequence(raw_seq)

            reason = failing_reason(scored)
            if reason is not None:
                exclusion_counts[reason] = exclusion_counts.get(reason, 0) + 1
                continue

            out_row = {"id": seq_id, "sequence": clean_sequence(raw_seq)}
            out_row.update(scored)

            if writer is None:
                writer = csv.DictWriter(out_f, fieldnames=list(out_row.keys()))
                writer.writeheader()
            writer.writerow(out_row)
            n_passed += 1

    out_f.close()

    print()
    print(f"Total rows read:     {n_total:,}")
    print(f"Empty sequence rows: {n_empty:,}")
    print(f"Passed all filters:  {n_passed:,}  -> {out_path}")
    print("Excluded breakdown:")
    for reason, count in sorted(exclusion_counts.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count:,}")


if __name__ == "__main__":
    main()