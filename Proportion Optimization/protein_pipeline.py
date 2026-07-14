import math
import random
import os
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

try:
    from localcider.sequenceParameters import SequenceParameters
    HAS_LOCALCIDER = True
except Exception:
    SequenceParameters = None
    HAS_LOCALCIDER = False

# Cache for SequenceParameters to avoid recreation overhead
# Using a simple dict with size limit to prevent memory bloat
CIDER_CACHE = {}
CIDER_CACHE_MAX_SIZE = 10000  # Adjust based on memory constraints

PKA = {"D": 3.9, "E": 4.1, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5, "Nterm": 8.0, "Cterm": 3.1}
VALID_AA = "ACDEFGHIKLMNPQRSTVWY"
CORE_ANALYSIS_COLUMNS = [
    "length",
    "molecular_weight",
    "isoelectric_point",
    "gravy",
    "aromaticity",
    "instability_index",
    "stable",
    "charge_at_pH7",
    "positive_res_RK",
    "negative_res_DE",
    "charge_sum",
    "helix_fraction",
    "sheet_fraction",
    "turn_fraction",
    "helix_sheet_sum",
]
EXTENDED_CHARGE_COLUMNS = ["fcr", "ncpr"]
CIDER_ANALYSIS_COLUMNS = ["cider_kappa"]
STATUS_COLUMNS = ["error"]

# ── CAHS motif-1 helix/amphipathicity scoring ────────────────────────────────
# Ported from the CAHS motif1 filter pipeline. Scores a fixed helical motif
# (construct residue numbering below) for helix propensity, hydrophobic
# moment / amphipathic face occupancy, salt-bridge count, and Pro/Gly
# "helix-breaker" content.
CORE_START_RESNUM = 121
MOTIF1_START_RESNUM = 124
MOTIF1_END_RESNUM = 142

# Pace & Scholtz (1998) helix propensity scale, kcal/mol (lower = more helix-favorable).
PACE_SCHOLTZ = {
    "A": 0.00, "R": 0.21, "N": 0.65, "D": 0.69, "C": 0.68,
    "Q": 0.39, "E": 0.40, "G": 1.00, "H": 0.61, "I": 0.41,
    "L": 0.21, "K": 0.26, "M": 0.24, "F": 0.54, "P": 3.16,
    "S": 0.50, "T": 0.66, "W": 0.49, "Y": 0.53, "V": 0.61,
}

# Eisenberg consensus hydrophobicity scale, used for the helical hydrophobic moment.
EISENBERG = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

CAHS_MOTIF1_DELTA_DEG = 100.0  # alpha-helix residues-per-turn angle

CAHS_MOTIF1_COLUMNS = [
    "motif1_sequence",
    "pace_scholtz_sum",
    "pace_scholtz_mean",
    "hydrophobic_moment",
    "helical_face_occupancy",
    "salt_bridge_count",
    "num_proline",
    "num_glycine",
    "helix_breaker_flag",
]

# ── Conservative substitution groups ─────────────────────────────────────────
CONSERVATIVE_GROUPS = {
    "G": ["A"],
    "A": ["G", "V", "S"],
    "V": ["A", "I", "L"],
    "L": ["V", "I", "M"],
    "I": ["V", "L", "M"],
    "M": ["L", "I"],
    "F": ["Y", "W"],
    "Y": ["F", "W", "H"],
    "W": ["F", "Y"],
    "S": ["T", "A", "N"],
    "T": ["S", "V", "N"],
    "N": ["Q", "S", "D"],
    "Q": ["N", "E", "K"],
    "K": ["R", "Q"],
    "R": ["K", "H"],
    "H": ["R", "K", "Y"],
    "D": ["E", "N"],
    "E": ["D", "Q"],
    "C": ["S", "A"],
    "P": ["A", "G"],
}

def extract_motif1(full_seq, motif_start=MOTIF1_START_RESNUM, motif_end=MOTIF1_END_RESNUM,
                    core_start=CORE_START_RESNUM):
    """Slice out the motif-1 helical segment using construct residue numbering."""
    start_idx = motif_start - core_start
    end_idx = motif_end - core_start + 1
    return full_seq[start_idx:end_idx]


def pace_scholtz_scores(seq):
    """Sum and mean Pace & Scholtz helix-propensity score across a sequence."""
    vals = [PACE_SCHOLTZ.get(aa) for aa in seq]
    valid = [v for v in vals if v is not None]
    if not valid:
        return None, None
    return sum(valid), sum(valid) / len(valid)


def hydrophobic_moment_and_face(seq, delta_deg=CAHS_MOTIF1_DELTA_DEG):
    """
    Eisenberg hydrophobic moment (normalized by length) and the fraction of
    hydrophobic residues sitting within 60 degrees of the moment vector
    ("helical face occupancy").
    """
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


def salt_bridge_score(seq):
    """Count i,i+3 / i,i+4 acidic-basic pairs (classic helical salt bridges)."""
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


def helix_breaker_flags(seq):
    """Pro/Gly counts and a combined hard-filter flag (either one breaks helix)."""
    n_pro, n_gly = seq.count("P"), seq.count("G")
    return n_pro, n_gly, (n_pro > 0 or n_gly > 0)


def analyze_cahs_motif1(full_seq, motif_start=MOTIF1_START_RESNUM, motif_end=MOTIF1_END_RESNUM,
                         core_start=CORE_START_RESNUM):
    """
    Score the CAHS motif-1 helical segment of a full-length sequence.
    Returns a dict matching CAHS_MOTIF1_COLUMNS. `full_seq` should already be
    cleaned (uppercase, no whitespace) -- analyze_sequence() does this before
    calling in.
    """
    motif1_seq = extract_motif1(full_seq, motif_start=motif_start, motif_end=motif_end,
                                 core_start=core_start)

    ps_sum, ps_mean = pace_scholtz_scores(motif1_seq)
    moment, face_score = hydrophobic_moment_and_face(motif1_seq)
    sb_score = salt_bridge_score(motif1_seq)
    n_pro, n_gly, breaker_flag = helix_breaker_flags(motif1_seq)

    return {
        "motif1_sequence": motif1_seq,
        "pace_scholtz_sum": round(ps_sum, 3) if ps_sum is not None else None,
        "pace_scholtz_mean": round(ps_mean, 3) if ps_mean is not None else None,
        "hydrophobic_moment": round(moment, 3) if moment is not None else None,
        "helical_face_occupancy": round(face_score, 3) if face_score is not None else None,
        "salt_bridge_count": sb_score,
        "num_proline": n_pro,
        "num_glycine": n_gly,
        "helix_breaker_flag": breaker_flag,
    }


def read_csv_safe(path_or_buffer, **kwargs):
    """
    Read a CSV while silently stripping trailing commas/whitespace that Excel
    and Google Sheets commonly append.  Drop any columns that are entirely
    empty after stripping (Unnamed: N artefacts).
    """
    import io
    if hasattr(path_or_buffer, "read"):          # file-like object (Streamlit upload)
        raw = path_or_buffer.read()
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8-sig")
    else:                                         # file path string
        with open(path_or_buffer, newline="", encoding="utf-8-sig") as f:
            raw = f.read()

    cleaned = "\n".join(line.rstrip().rstrip(",") for line in raw.splitlines())
    df = pd.read_csv(io.StringIO(cleaned), **kwargs)

    # Drop any phantom columns from leftover commas in the header itself
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    if unnamed:
        df = df.drop(columns=unnamed)
    all_null = [c for c in df.columns if df[c].isna().all()]
    if all_null:
        df = df.drop(columns=all_null)
    return df
    # Drop any phantom columns from leftover commas in the header itself
    unnamed = [c for c in df.columns if str(c).startswith("Unnamed:")]
    if unnamed:
        df = df.drop(columns=unnamed)
    all_null = [c for c in df.columns if df[c].isna().all()]
    if all_null:
        df = df.drop(columns=all_null)
    return df
def net_charge(seq, ph=7.0, counts=None):
    if counts is None:
        counts = Counter(seq)
    c = 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"]))
    c -= 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph))
    for aa in ("D", "E", "C", "Y"):
        c -= counts[aa] / (1.0 + 10 ** (PKA[aa] - ph))
    for aa in ("H", "K", "R"):
        c += counts[aa] / (1.0 + 10 ** (ph - PKA[aa]))
    return round(c, 3)


def calculate_fcr(seq, counts=None):
    """
    Fractional Charged Residues.
    FCR = (K + R + E + D + H_ionized) / length
    H is conditionally included (ionized at pH 7): ~1 - 1/(1 + 10^(pH - pKa_H))
    At pH 7 with pKa=6.0, H is ~10% ionized, so weight as 0.1
    """
    n = len(seq)
    if n == 0:
        return 0.0
    if counts is None:
        counts = Counter(seq)
    positive = counts["K"] + counts["R"]
    negative = counts["E"] + counts["D"]
    h_fraction = 0.1  # H ionization at pH 7 (pKa 6.0)
    h_contribution = counts["H"] * h_fraction
    fcr = (positive + negative + h_contribution) / n
    return round(fcr, 4)


def calculate_ncpr(seq, counts=None):
    """
    Net Charge Per Residue.
    NCPR = (K + R + H_ionized - E - D) / length
    H contributes +0.1 at pH 7 (≈10% protonated)
    """
    n = len(seq)
    if n == 0:
        return 0.0
    if counts is None:
        counts = Counter(seq)
    positive = counts["K"] + counts["R"]
    negative = counts["E"] + counts["D"]
    h_fraction = 0.1
    h_contribution = counts["H"] * h_fraction
    ncpr = (positive + h_contribution - negative) / n
    return round(ncpr, 4)


def calculate_cider_kappa(seq):
    """
    CIDER κ (kappa) charge-patterning metric from Das & Pappu (2013).
    Requires localcider: pip install localcider
    Uses caching to avoid recreating SequenceParameters for repeated sequences.
    """
    if not seq or not HAS_LOCALCIDER:
        return None

    # Check cache first
    if seq in CIDER_CACHE:
        sp = CIDER_CACHE[seq]
    else:
        try:
            sp = SequenceParameters(seq)
            # Simple cache management - if full, clear it (could be more sophisticated)
            if len(CIDER_CACHE) >= CIDER_CACHE_MAX_SIZE:
                CIDER_CACHE.clear()
            CIDER_CACHE[seq] = sp
        except Exception:
            return None

    try:
        return round(sp.get_kappa(), 4)
    except Exception:
        return None

def analyze_sequence(seq, include_extended_charge_metrics=True, include_cider_kappa=True,
                      include_cahs_motif1=False):
    seq = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(seq) - set(VALID_AA)
    if bad:
        return {"error": f"invalid chars: {bad}"}
    if len(seq) < 2:
        return {"error": "too short"}
    try:
        counts = Counter(seq)
        pa = ProteinAnalysis(seq)
        helix, turn, sheet = pa.secondary_structure_fraction()
        positive_res = counts["R"] + counts["K"]
        negative_res = counts["D"] + counts["E"]
        helix_sheet_sum = helix + sheet
        result = {
            "length": len(seq),
            "molecular_weight": round(pa.molecular_weight(), 2),
            "isoelectric_point": round(pa.isoelectric_point(), 2),
            "gravy": round(pa.gravy(), 4),
            "aromaticity": round(pa.aromaticity(), 4),
            "instability_index": round(pa.instability_index(), 2),
            "stable": pa.instability_index() <= 40,
            "charge_at_pH7": net_charge(seq, counts=counts),
            "positive_res_RK": positive_res,
            "negative_res_DE": negative_res,
            "charge_sum": positive_res + negative_res,
            "helix_fraction": round(helix, 4),
            "sheet_fraction": round(sheet, 4),
            "turn_fraction": round(turn, 4),
            "helix_sheet_sum": round(helix_sheet_sum, 4),
            "fcr": calculate_fcr(seq, counts=counts) if include_extended_charge_metrics else None,
            "ncpr": calculate_ncpr(seq, counts=counts) if include_extended_charge_metrics else None,
            "cider_kappa": calculate_cider_kappa(seq) if include_cider_kappa else None,
            "error": "",
        }
        if include_cahs_motif1:
            try:
                result.update(analyze_cahs_motif1(seq))
            except Exception as e:
                result["error"] = f"cahs_motif1: {e}"
        return result
    except Exception as e:
        return {"error": str(e)}


def get_requested_analysis_columns(include_extended_charge_metrics=True, include_cider_kappa=True,
                                    include_cahs_motif1=False):
    columns = list(CORE_ANALYSIS_COLUMNS)
    if include_extended_charge_metrics:
        columns.extend(EXTENDED_CHARGE_COLUMNS)
    if include_cider_kappa:
        columns.extend(CIDER_ANALYSIS_COLUMNS)
    if include_cahs_motif1:
        columns.extend(CAHS_MOTIF1_COLUMNS)
    columns.extend(STATUS_COLUMNS)
    return columns


def _empty_value_mask(series, treat_blank_string_as_empty=True):
    mask = series.isna()
    if pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series):
        text = series.astype("string")
        normalized = text.str.strip().str.lower()
        empty_tokens = ["nan", "none"]
        if treat_blank_string_as_empty:
            empty_tokens.insert(0, "")
        mask = mask | normalized.isin(empty_tokens)
    return mask.fillna(True)


def screen_analysis_columns(
    df,
    include_extended_charge_metrics=True,
    include_cider_kappa=True,
    include_cahs_motif1=False,
    force_recompute_mask=None,
):
    """
    Inspect which analyzer output columns are already available and which rows
    still need to be computed.

    Returns a dict with requested columns, missing columns, partially empty
    columns, and a row mask indicating which rows need fresh analysis.
    """
    requested_columns = get_requested_analysis_columns(
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )

    if force_recompute_mask is None:
        force_recompute_mask = pd.Series(False, index=df.index)
    else:
        force_recompute_mask = pd.Series(force_recompute_mask, index=df.index).fillna(False).astype(bool)

    rows_to_update_mask = force_recompute_mask.copy()
    missing_columns = []
    partially_empty_columns = []
    empty_value_counts = {}

    for col in requested_columns:
        should_trigger_update = col not in STATUS_COLUMNS

        if col not in df.columns:
            if should_trigger_update:
                missing_columns.append(col)
            empty_mask = pd.Series(True, index=df.index)
            if should_trigger_update:
                empty_value_counts[col] = len(df)
        else:
            empty_mask = _empty_value_mask(df[col], treat_blank_string_as_empty=(col not in STATUS_COLUMNS))
            empty_count = int(empty_mask.sum())
            if empty_count and should_trigger_update:
                partially_empty_columns.append(col)
                empty_value_counts[col] = empty_count

        if should_trigger_update:
            rows_to_update_mask |= empty_mask

    return {
        "requested_columns": requested_columns,
        "missing_columns": missing_columns,
        "partially_empty_columns": partially_empty_columns,
        "empty_value_counts": empty_value_counts,
        "rows_to_update_mask": rows_to_update_mask,
        "rows_to_update_count": int(rows_to_update_mask.sum()),
        "forced_recompute_count": int(force_recompute_mask.sum()),
    }


def ensure_analysis_columns(
    df,
    include_extended_charge_metrics=True,
    include_cider_kappa=True,
    include_cahs_motif1=False,
):
    """
    Ensure requested analyzer columns exist in a dataframe. Missing status
    columns are initialized to blank strings; missing metric columns are
    initialized to pandas NA.
    """
    requested_columns = get_requested_analysis_columns(
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )

    normalized_df = df.copy()
    for col in requested_columns:
        if col not in normalized_df.columns:
            normalized_df[col] = "" if col in STATUS_COLUMNS else pd.NA

    for col in STATUS_COLUMNS:
        if col in normalized_df.columns:
            normalized_df[col] = normalized_df[col].astype(object)
            normalized_df.loc[normalized_df[col].isna(), col] = ""

    return normalized_df


def merge_analysis_results(
    df,
    results_df,
    include_extended_charge_metrics=True,
    include_cider_kappa=True,
    include_cahs_motif1=False,
    force_recompute_mask=None,
):
    """
    Merge fresh analyzer results into an existing dataframe, filling blank
    analysis cells and adding any missing analyzer columns without duplicating
    them.
    """
    requested_columns = get_requested_analysis_columns(
        include_extended_charge_metrics=include_extended_charge_metrics,
        include_cider_kappa=include_cider_kappa,
        include_cahs_motif1=include_cahs_motif1,
    )

    if force_recompute_mask is None:
        force_recompute_mask = pd.Series(False, index=df.index)
    else:
        force_recompute_mask = pd.Series(force_recompute_mask, index=df.index).fillna(False).astype(bool)

    merged_df = df.copy()

    for col in requested_columns:
        if col not in merged_df.columns:
            merged_df[col] = pd.NA
        if col not in results_df.columns:
            continue

        if col in STATUS_COLUMNS:
            merged_df[col] = merged_df[col].astype(object)

        incoming = results_df[col].reindex(merged_df.index)
        fill_mask = _empty_value_mask(
            merged_df[col],
            treat_blank_string_as_empty=(col not in STATUS_COLUMNS),
        ) | force_recompute_mask
        fill_mask &= incoming.notna()
        merged_df.loc[fill_mask, col] = incoming.loc[fill_mask]

    return merged_df


def parse_mutation_regions(text):
    regions = []
    if not str(text).strip():
        return regions
    parts = [p.strip() for p in str(text).split(",") if p.strip()]
    for p in parts:
        if "-" in p:
            a, b = [x.strip() for x in p.split("-", 1)]
            if not a.isdigit() or not b.isdigit():
                raise ValueError(f"Invalid region: {p}")
            start, end = int(a), int(b)
        else:
            if not p.isdigit():
                raise ValueError(f"Invalid position: {p}")
            start = end = int(p)
        if start < 1 or end < 1:
            raise ValueError(f"Positions must be >= 1: {p}")
        if end < start:
            raise ValueError(f"Region end must be >= start: {p}")
        regions.append((start, end))
    return regions


def parse_copy_count(value):
    if pd.isna(value) or str(value).strip() == "":
        return 1
    s = str(value).strip()
    try:
        num = float(s)
    except ValueError as e:
        raise ValueError(f"Invalid copy count: {s}") from e
    if not num.is_integer():
        raise ValueError(f"Copy count must be an integer: {s}")
    n = int(num)
    if n < 0:
        raise ValueError(f"Copy count must be >= 0: {s}")
    return n


# ── Combo count helpers ───────────────────────────────────────────────────────
def max_unique_conservative_variants(seq, regions):
    """
    Returns the maximum number of unique sequences producible by conservative
    substitution across the given regions (product of substitute pool sizes,
    NOT including the original residue — since we always mutate away from it).
    Used to detect when num_copies > unique possibilities.
    Returns (max_unique, per_position) where per_position is a list of
    (pos, original_aa, n_substitutes).
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    total = 1
    per_position = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            aa = clean[idx]
            subs = CONSERVATIVE_GROUPS.get(aa, [])
            n = len(subs) if subs else 1  # if no subs, position is fixed → factor of 1
            per_position.append((pos, aa, n))
            total *= n
    return total, per_position


def count_conservative_combos(seq, regions):
    """
    Returns (total_combos, per_position_breakdown).
    total_combos = product of (n_substitutes + 1) across all positions in region,
    where +1 counts keeping the original residue as an option.
    per_position_breakdown = list of (position, original_aa, substitutes, n_choices).
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    total = 1
    per_position = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            aa = clean[idx]
            subs = CONSERVATIVE_GROUPS.get(aa, [])
            n_choices = len(subs) + 1  # include original
            per_position.append((pos, aa, subs, n_choices))
            total *= n_choices
    return total, per_position


def count_random_combos(seq, regions):
    """
    Returns (total_combos, n_positions).
    Each position can be any of 20 AAs, so total = 20^n_positions.
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    seen = set()
    n_positions = 0
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if 0 <= idx < len(clean):
                n_positions += 1
    return 20 ** n_positions, n_positions


# ── Mode 1: random substitution ───────────────────────────────────────────────
def mutate_sequence_random(seq, regions):
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")
    arr = list(clean)
    changed_positions = []
    for start, end in regions:
        for pos in range(start, end + 1):
            idx = pos - 1
            if idx < 0 or idx >= len(arr):
                continue
            current = arr[idx]
            choices = [aa for aa in VALID_AA if aa != current]
            arr[idx] = random.choice(choices)
            changed_positions.append(pos)
    return "".join(arr), len(changed_positions)


# ── Mode 2: conservative substitution ────────────────────────────────────────
def mutate_sequence_conservative(seq, regions):
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")
    arr = list(clean)
    changed_positions = []
    skipped_positions = []
    for start, end in regions:
        for pos in range(start, end + 1):
            idx = pos - 1
            if idx < 0 or idx >= len(arr):
                continue
            current = arr[idx]
            choices = CONSERVATIVE_GROUPS.get(current, [])
            if choices:
                arr[idx] = random.choice(choices)
                changed_positions.append(pos)
            else:
                skipped_positions.append(pos)
    return "".join(arr), len(changed_positions), skipped_positions


# ── Mode 3: exhaustive single-position scanning ───────────────────────────────
def scan_sequence_single_position(seq, regions):
    """
    For each position in regions, generate one mutant per conservative substitute,
    mutating only that single position at a time.
    Returns list of {mutated_seq, mutated_position, original_aa, new_aa}.
    """
    clean = str(seq).strip().upper().replace(" ", "").replace("\n", "")
    bad = set(clean) - set(VALID_AA)
    if bad:
        raise ValueError(f"invalid chars: {bad}")
    if not clean:
        raise ValueError("empty sequence")

    variants = []
    seen = set()
    for start, end in regions:
        for pos in range(start, end + 1):
            if pos in seen:
                continue
            seen.add(pos)
            idx = pos - 1
            if idx < 0 or idx >= len(clean):
                continue
            original_aa = clean[idx]
            for new_aa in CONSERVATIVE_GROUPS.get(original_aa, []):
                arr = list(clean)
                arr[idx] = new_aa
                variants.append({
                    "mutated_seq": "".join(arr),
                    "mutated_position": pos,
                    "original_aa": original_aa,
                    "new_aa": new_aa,
                })
    return variants


# ── Shared expansion logic ────────────────────────────────────────────────────
def expand_rows_with_mutations(df, seq_col, region_col, copies_col,
                               random_seed, mutation_mode="random"):
    """
    mutation_mode: "random" | "conservative" | "scan"
    In scan mode, copies_col is ignored — all single-position variants are
    generated exhaustively.

    Returns (expanded_df, skipped_zero_copy_rows, duplicate_warnings)
    where duplicate_warnings is a list of dicts:
        {row: int, requested: int, max_unique: int, per_position: list}
    """
    random.seed(int(random_seed))
    expanded_rows = []
    skipped_zero_copy_rows = 0
    duplicate_warnings = []

    for idx, row in df.iterrows():
        raw_seq = row[seq_col]
        raw_region = row[region_col]
        raw_copies = row[copies_col]

        region_text = "" if pd.isna(raw_region) else str(raw_region).strip()
        row_error = ""

        try:
            regions = parse_mutation_regions(region_text)
        except ValueError as e:
            regions = []
            row_error = str(e)

        # ── Scan mode ─────────────────────────────────────────────────────────
        if mutation_mode == "scan":
            base_row = row.copy()
            base_row["source_row"] = idx + 1
            base_row["original_sequence"] = str(raw_seq)
            base_row["mutation_regions"] = region_text
            base_row["mutation_mode"] = mutation_mode

            if row_error:
                base_row["mutation_error"] = row_error
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            try:
                variants = scan_sequence_single_position(raw_seq, regions)
            except ValueError as e:
                base_row["mutation_error"] = str(e)
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            if not variants:
                base_row["mutation_error"] = "no conservative substitutes found in region"
                base_row["mutated_position"] = ""
                base_row["original_aa"] = ""
                base_row["new_aa"] = ""
                expanded_rows.append(base_row)
                continue

            for v in variants:
                new_row = base_row.copy()
                new_row[seq_col] = v["mutated_seq"]
                new_row["mutated_position"] = v["mutated_position"]
                new_row["original_aa"] = v["original_aa"]
                new_row["new_aa"] = v["new_aa"]
                new_row["mutation_error"] = ""
                expanded_rows.append(new_row)
            continue

        # ── Random / conservative modes ───────────────────────────────────────
        try:
            num_copies = parse_copy_count(raw_copies)
        except ValueError as e:
            num_copies = 1
            row_error = str(e) if not row_error else f"{row_error}; {e}"

        if num_copies == 0:
            skipped_zero_copy_rows += 1
            continue

        # ── Duplicate check for conservative and random modes ─────────────────
        if not row_error and regions:
            try:
                if mutation_mode == "conservative":
                    max_unique, per_pos = max_unique_conservative_variants(str(raw_seq), regions)
                    if num_copies > max_unique:
                        duplicate_warnings.append({
                            "row": idx + 1,
                            "requested": num_copies,
                            "max_unique": max_unique,
                            "mode": "conservative",
                            "per_position": per_pos,
                        })
                elif mutation_mode == "random":
                    max_unique, n_positions = count_random_combos(str(raw_seq), regions)
                    if num_copies > max_unique:
                        duplicate_warnings.append({
                            "row": idx + 1,
                            "requested": num_copies,
                            "max_unique": max_unique,
                            "mode": "random",
                            "per_position": [(None, None, 19)] * n_positions,
                        })
            except Exception:
                pass

        for copy_i in range(1, num_copies + 1):
            new_row = row.copy()
            new_row["source_row"] = idx + 1
            new_row["copy_index"] = copy_i
            new_row["requested_copies"] = num_copies
            new_row["original_sequence"] = str(raw_seq)
            new_row["mutation_regions"] = region_text
            new_row["mutation_mode"] = mutation_mode

            if row_error:
                new_row[seq_col] = str(raw_seq)
                new_row["mutated_positions_count"] = 0
                new_row["conservative_skipped"] = ""
                new_row["mutation_error"] = row_error
            else:
                try:
                    if mutation_mode == "conservative":
                        mutated_seq, n_changed, skipped = mutate_sequence_conservative(raw_seq, regions)
                        new_row["conservative_skipped"] = ",".join(map(str, skipped)) if skipped else ""
                    else:
                        mutated_seq, n_changed = mutate_sequence_random(raw_seq, regions)
                        new_row["conservative_skipped"] = ""
                    new_row[seq_col] = mutated_seq
                    new_row["mutated_positions_count"] = n_changed
                    new_row["mutation_error"] = ""
                except ValueError as e:
                    new_row[seq_col] = str(raw_seq)
                    new_row["mutated_positions_count"] = 0
                    new_row["conservative_skipped"] = ""
                    new_row["mutation_error"] = str(e)

            expanded_rows.append(new_row)

    return pd.DataFrame(expanded_rows), skipped_zero_copy_rows, duplicate_warnings


# ── Parallel analysis ────────────────────────────────────────────────────────
def analyze_sequences_parallel(
    sequences,
    num_workers=None,
    progress_callback=None,
    use_processes=False,
    include_extended_charge_metrics=True,
    include_cider_kappa=True,
    include_cahs_motif1=False,
    max_pending_tasks=None,
):
    """
    Analyze multiple sequences in parallel.

    Args:
        sequences: List or pandas Series of sequence strings
        num_workers: Number of worker threads/processes (default: CPU count)
        progress_callback: Optional callable(completed, total) for progress updates
        use_processes: Use ProcessPoolExecutor when True; ThreadPoolExecutor when False (default)
        include_extended_charge_metrics: Include FCR and NCPR in output
        include_cider_kappa: Include localCIDER kappa in output (slowest metric)
        include_cahs_motif1: Include CAHS motif-1 helix/amphipathicity scoring in output
        max_pending_tasks: Max number of in-flight futures to keep queued

    Returns:
        List of analysis dictionaries
    """
    if num_workers is None:
        num_workers = (os.cpu_count() or 4) if use_processes else min(os.cpu_count() or 4, 8)
    if num_workers < 1:
        num_workers = 1

    total = len(sequences)
    if total == 0:
        if progress_callback:
            progress_callback(0, 0)
        return []

    if max_pending_tasks is None:
        pending_factor = 2 if use_processes else 8
        max_pending_tasks = max(num_workers * pending_factor, num_workers)
    else:
        max_pending_tasks = max(1, int(max_pending_tasks))

    executor_cls = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    results = [None] * total
    completed = 0

    with executor_cls(max_workers=num_workers) as executor:
        future_to_index = {}
        sequence_iter = iter(enumerate(sequences))

        def submit_next():
            try:
                idx, seq = next(sequence_iter)
            except StopIteration:
                return False

            future = executor.submit(
                analyze_sequence,
                str(seq),
                include_extended_charge_metrics=include_extended_charge_metrics,
                include_cider_kappa=include_cider_kappa,
                include_cahs_motif1=include_cahs_motif1,
            )
            future_to_index[future] = idx
            return True

        initial_submissions = min(total, max_pending_tasks)
        for _ in range(initial_submissions):
            submit_next()

        while future_to_index:
            done, _ = wait(future_to_index, return_when=FIRST_COMPLETED)
            for future in done:
                idx = future_to_index.pop(future)
                results[idx] = future.result()
                completed += 1
                if progress_callback:
                    progress_callback(completed, total)
                submit_next()

    return results