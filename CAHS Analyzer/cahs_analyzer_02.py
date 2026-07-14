#!/usr/bin/env python3
"""
motif1_helix_scorer.py

Batch-scores a CSV of protein/peptide sequences (e.g. CAHS motif 1 mutants)
on the parameter set discussed for filtering a directed-evolution library:

  Helical propensity   -> Pace & Scholtz (1998) ddG helix-propensity scale
  Amphipathicity       -> Eisenberg hydrophobic moment
  Helical-face check   -> fraction of hydrophobic residues clustered near
                          the hydrophobic-moment vector (helical-wheel face)
  Salt-bridge spacing  -> count of i,i+3 / i,i+4 E/K and D/R ion pairs
  Hard filters         -> Pro/Gly count (helix-breaking flags)
  Baseline composition -> GRAVY, instability index, aliphatic index,
                          net charge at pH 7, isoelectric point
                          (via Biopython ProtParam, same numbers your
                          existing Streamlit tool reports)
  Optional             -> BLOSUM62 similarity score vs a wild-type
                          reference sequence, if you set WILDTYPE below

NOTE ON DATA PROVENANCE:
- Pace & Scholtz values are taken from Pace & Scholtz, Biophys J. 1998;
  75(1):422-7 (verified against the published table).
- Eisenberg consensus hydrophobicity values are the commonly reproduced
  Eisenberg (1984) consensus scale. Double check both tables against the
  original papers before using results in anything you publish -- these
  are screening heuristics, not a replacement for AGADIR / experimental
  CD / IUPred3, which you should still run on your final shortlist.

USAGE:
    python motif1_helix_scorer.py input.csv output.csv

INPUT CSV:
    Must have a column named "sequence" (case-insensitive).
    An optional "id" column will be carried through; if absent, sequences
    are numbered seq_1, seq_2, ...

OUTPUT CSV:
    One row per input sequence with all scored columns appended.
"""

import sys
import csv
import math

# ---------------------------------------------------------------------------
# CONSTRUCT NUMBERING
# Set these to match your CAHS Core construct. CORE_START_RESNUM is the
# residue number corresponding to the FIRST character of every sequence in
# your input CSV. MOTIF1_START/END_RESNUM are the (inclusive) residue
# numbers of motif 1 within that same numbering scheme. The script uses
# these to automatically slice out motif 1 from each full-length row.
#
# Example based on your CAHS Core (121-183):
#   CORE_START_RESNUM   = 121
#   MOTIF1_START_RESNUM = 124
#   MOTIF1_END_RESNUM   = 142
# ---------------------------------------------------------------------------
CORE_START_RESNUM = 121
MOTIF1_START_RESNUM = 124
MOTIF1_END_RESNUM = 142

# ---------------------------------------------------------------------------
# Optional: set the FULL wild-type core sequence here (same length/frame as
# your mutant rows) to get, per mutant: a whole-sequence BLOSUM62 similarity
# score, plus deltas of the motif-1-only helix propensity and hydrophobic
# moment relative to wild-type motif 1. Leave as None to skip these columns.
# ---------------------------------------------------------------------------
FULL_WILDTYPE = None  # e.g. "TEAYRKQQEVEADKIRKELEKQHLRDVEFRKDIVEMAIENQKKMIDVESRYAKKDMDRERVKV"

# Thresholds for the two-sided review flags (tune these after looking at the
# actual spread of scores in your own batch -- these are starting points,
# not literature-derived cutoffs).
OVERRIGIDITY_PS_DELTA_THRESHOLD = -0.15   # mutant motif1 mean more helix-favorable than WT by this much (kcal/mol) -> flag for review
EXTREME_MOMENT_RATIO_THRESHOLD = 1.5      # mutant moment > this multiple of WT moment -> flag for review

# ---------------------------------------------------------------------------
# Pace & Scholtz (1998) helix propensity scale, kcal/mol (LOWER = more
# helix-favorable; Ala = 0 is the reference, Pro/Gly are the worst).
# ---------------------------------------------------------------------------
PACE_SCHOLTZ = {
    "A": 0.00, "R": 0.21, "N": 0.65, "D": 0.69, "C": 0.68,
    "Q": 0.39, "E": 0.40, "G": 1.00, "H": 0.61, "I": 0.41,
    "L": 0.21, "K": 0.26, "M": 0.24, "F": 0.54, "P": 3.16,
    "S": 0.50, "T": 0.66, "W": 0.49, "Y": 0.53, "V": 0.61,
}

# ---------------------------------------------------------------------------
# Eisenberg consensus hydrophobicity scale (used for hydrophobic moment).
# ---------------------------------------------------------------------------
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

DELTA_DEG = 100.0  # degrees per residue for an alpha helix (3.6 res/turn)

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.Align import substitution_matrices
    BLOSUM62 = substitution_matrices.load("BLOSUM62")
    HAVE_BIOPYTHON = True
except ImportError:
    HAVE_BIOPYTHON = False


def clean_sequence(seq: str) -> str:
    return "".join(seq.upper().split())


def pace_scholtz_scores(seq: str):
    vals = [PACE_SCHOLTZ.get(aa) for aa in seq]
    valid = [v for v in vals if v is not None]
    if not valid:
        return None, None
    return sum(valid), sum(valid) / len(valid)


def hydrophobic_moment_and_face(seq: str, delta_deg: float = DELTA_DEG):
    """
    Eisenberg hydrophobic moment (magnitude, normalized by length) plus a
    'face occupancy' score: the fraction of hydrophobic residues (Eisenberg
    value > 0) whose angular position on the helical wheel falls within
    +/-60 degrees of the net hydrophobic moment vector. High magnitude +
    high face occupancy = a clean, well-registered amphipathic helix.
    Low face occupancy despite high magnitude flags a mutation that raised
    average hydrophobicity without keeping residues on the same face.
    """
    delta_rad = math.radians(delta_deg)
    sin_sum, cos_sum = 0.0, 0.0
    angles = []
    hydrophobic_flags = []

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

    on_face = 0
    total_hydrophobic = 0
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


def salt_bridge_score(seq: str):
    """
    Counts i,i+3 and i,i+4 E/K and D/R ion-pair spacings in either
    direction (the classic i,i+3/i,i+4 helix-stabilizing salt bridge
    geometry). Returns a raw count -- higher is generally more
    helix-stabilizing, but check spacing makes biological sense in
    context (this is a cheap proxy, not a real electrostatics calc).
    """
    acidic = {"E", "D"}
    basic = {"K", "R"}
    count = 0
    n = len(seq)
    for i in range(n):
        for offset in (3, 4):
            j = i + offset
            if j >= n:
                continue
            pair = {seq[i], seq[j]}
            if (seq[i] in acidic and seq[j] in basic) or (seq[i] in basic and seq[j] in acidic):
                count += 1
    return count


def helix_breaker_flags(seq: str):
    n_pro = seq.count("P")
    n_gly = seq.count("G")
    return n_pro, n_gly, (n_pro > 0 or n_gly > 0)


def blosum62_score(seq: str, ref: str):
    if not ref or len(seq) != len(ref):
        return None
    score = 0
    for a, b in zip(seq, ref):
        try:
            score += BLOSUM62[a, b]
        except KeyError:
            return None
    return score


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


def extract_motif1(full_seq: str):
    start_idx = MOTIF1_START_RESNUM - CORE_START_RESNUM
    end_idx = MOTIF1_END_RESNUM - CORE_START_RESNUM + 1  # +1: slice end is exclusive
    # print(full_seq[start_idx:end_idx])
    return full_seq[start_idx:end_idx]


# Pre-compute wild-type motif 1 reference metrics once, if provided.
_WT_MOTIF1 = None
_WT_PS_MEAN = None
_WT_MOMENT = None
if FULL_WILDTYPE:
    _WT_MOTIF1 = extract_motif1(clean_sequence(FULL_WILDTYPE))
    _, _WT_PS_MEAN = pace_scholtz_scores(_WT_MOTIF1)
    _WT_MOMENT, _ = hydrophobic_moment_and_face(_WT_MOTIF1)


def score_sequence(full_seq: str):
    full_seq = clean_sequence(full_seq)
    motif1_seq = extract_motif1(full_seq)

    # --- motif-1-only helix-quality metrics ---
    ps_sum, ps_mean = pace_scholtz_scores(motif1_seq)
    moment, face_score = hydrophobic_moment_and_face(motif1_seq)
    sb_score = salt_bridge_score(motif1_seq)
    n_pro, n_gly, breaker_flag = helix_breaker_flags(motif1_seq)

    # --- whole-construct baseline (expression/solubility relevant) ---
    baseline = biopython_baseline(full_seq)

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
    row.update(baseline)

    if FULL_WILDTYPE:
        row["blosum62_vs_wildtype"] = blosum62_score(full_seq, clean_sequence(FULL_WILDTYPE))

        if ps_mean is not None and _WT_PS_MEAN is not None:
            ps_delta = ps_mean - _WT_PS_MEAN
            row["pace_scholtz_mean_delta_vs_wt"] = round(ps_delta, 3)
            row["flag_possible_overrigidity"] = ps_delta < OVERRIGIDITY_PS_DELTA_THRESHOLD
        else:
            row["pace_scholtz_mean_delta_vs_wt"] = None
            row["flag_possible_overrigidity"] = None

        if moment is not None and _WT_MOMENT:
            moment_ratio = moment / _WT_MOMENT
            row["hydrophobic_moment_ratio_vs_wt"] = round(moment_ratio, 3)
            row["flag_extreme_moment"] = moment_ratio > EXTREME_MOMENT_RATIO_THRESHOLD
        else:
            row["hydrophobic_moment_ratio_vs_wt"] = None
            row["flag_extreme_moment"] = None

    return row


def main():
    if len(sys.argv) != 3:
        print("Usage: python motif1_helix_scorer.py input.csv output.csv")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames_lower = {name.lower(): name for name in reader.fieldnames}
        if "sequence" not in fieldnames_lower:
            print("ERROR: input CSV must have a 'sequence' column.")
            sys.exit(1)
        seq_col = fieldnames_lower["sequence"]
        id_col = fieldnames_lower.get("id")

        rows_in = list(reader)

    results = []
    for i, row in enumerate(rows_in, start=1):
        seq_id = row.get(id_col, f"seq_{i}") if id_col else f"seq_{i}"
        raw_seq = row.get(seq_col, "")
        scored = score_sequence(raw_seq)
        out_row = {"id": seq_id, "sequence": clean_sequence(raw_seq)}
        out_row.update(scored)
        results.append(out_row)

    if not results:
        print("No sequences found in input CSV.")
        sys.exit(1)

    fieldnames = list(results[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Scored {len(results)} sequences -> {out_path}")


if __name__ == "__main__":
    main()