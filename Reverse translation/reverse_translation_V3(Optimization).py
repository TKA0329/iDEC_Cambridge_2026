"""
bio_model_v4.py
===============
Degenerate IUPAC codon design from a multiple sequence alignment (MSA).

Given an aligned FASTA file of protein sequences, the programme selects one
IUPAC degenerate codon per alignment column such that the set of amino acids
encoded by that codon best matches the amino acid diversity observed at that
position.  The result is a single degenerate DNA sequence that compactly
describes a protein library.

Three optimisation modes are available:

    direct   Directly minimises  2 − recall − precision  for each position.
             This is always the globally optimal solution for the stated
             objective and requires no user-tunable parameters.

    grid     Sweeps the two penalty weights (w_loss, w_offtarget) over a
             user-defined grid, evaluates the resulting  avg(2−R−P)  at each
             point, prints an ASCII heatmap of the landscape, and reports
             which weight pair best reproduces the direct optimum.

    surface  Runs the same grid, then fits a degree-2 polynomial response
             surface z = f(w_loss, w_offtarget) and finds its analytical
             minimum.  Reports whether that minimum agrees with direct mode
             and gives an R² quality-of-fit score.

Usage
-----
    python bio_model_v4.py <alignment.fasta> [options]

    --mode    direct | grid | surface   (default: direct)
    --grid-step F   step size for grid / surface sweep  (default: 0.1)
    --grid-max  F   upper bound of weight sweep          (default: 2.0)
    --verbose       print per-position detail

Examples
--------
    python bio_model_v4.py sequences.fasta
    python bio_model_v4.py sequences.fasta --mode grid --grid-step 0.2
    python bio_model_v4.py sequences.fasta --mode surface --verbose
"""

from collections import Counter
from itertools import product
import math
import argparse

# ═══════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════

GROUPS = {
    "hydrophobic": {"A","V","L","I","M","F","W","Y","P"},
    "polar":       {"S","T","N","Q","C","G"},
    "charged":     {"D","E","K","R","H"},
}

AA_TO_GROUP = {aa: g for g, aas in GROUPS.items() for aa in aas}

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# All 15 IUPAC nucleotide codes and the concrete bases they represent
IUPAC = {
    "A":{"A"},"C":{"C"},"G":{"G"},"T":{"T"},
    "R":{"A","G"},"Y":{"C","T"},"S":{"G","C"},
    "W":{"A","T"},"K":{"G","T"},"M":{"A","C"},
    "B":{"C","G","T"},"D":{"A","G","T"},
    "H":{"A","C","T"},"V":{"A","C","G"},"N":{"A","C","G","T"},
}


def _expand(codon):
    """Expand a 3-character IUPAC codon to all concrete DNA triplets."""
    return ("".join(t) for t in product(*[IUPAC[b] for b in codon]))


def _build_codon_aa_map():
    """
    Pre-compute {iupac_codon: frozenset_of_amino_acids} for all 15³ = 3375
    IUPAC codons.  Stop codons (*) are excluded.  Codons that encode only
    stop codons are dropped entirely.  Built once at import time.
    """
    mapping = {}
    for c in ("".join(t) for t in product(IUPAC.keys(), repeat=3)):
        aas = frozenset(
            CODON_TABLE[x]
            for x in _expand(c)
            if x in CODON_TABLE and CODON_TABLE[x] != "*"
        )
        if aas:
            mapping[c] = aas
    return mapping


# {iupac_codon: frozenset of amino acids}  — used by all search functions
CODON_AA_MAP = _build_codon_aa_map()


# ═══════════════════════════════════════════════════════════════
# I/O
# ═══════════════════════════════════════════════════════════════

def read_fasta(path):
    """
    Parse a FASTA file and return a list of uppercase sequence strings.
    Validates that all sequences have the same length (i.e. are aligned).
    """
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.upper())
    if cur:
        seqs.append("".join(cur))

    if not seqs:
        raise ValueError("No sequences found in FASTA file.")

    L = len(seqs[0])
    for i, s in enumerate(seqs[1:], start=2):
        if len(s) != L:
            raise ValueError(
                f"Sequence {i} has length {len(s)}, expected {L}. "
                "All sequences must be aligned to the same length."
            )
    return seqs


def columns(seqs):
    """
    Transpose aligned sequences into per-column lists.
    Gap characters ('-') are excluded from each column.
    """
    L = len(seqs[0])
    return [[s[i] for s in seqs if s[i] != "-"] for i in range(L)]


# ═══════════════════════════════════════════════════════════════
# Metrics
# ═══════════════════════════════════════════════════════════════

def shannon_entropy(col):
    """Shannon entropy of an amino acid column, in bits (log base 2)."""
    counts = Counter(col)
    total  = len(col)
    return -sum(
        (v / total) * math.log2(v / total)
        for v in counts.values() if v > 0
    )


def position_metrics(orig_set, enc_aas):
    """
    Compute all quality metrics for one alignment position.

    Parameters
    ----------
    orig_set : set   — amino acids observed in the alignment column
    enc_aas  : frozenset — amino acids encoded by the candidate IUPAC codon

    Returns  (as a dict)
    -------
    recall      : fraction of observed AAs that are encoded          (↑ better)
    precision   : fraction of encoded AAs that were observed         (↑ better)
    loss        : 1 − recall  — missed observed AAs                  (↓ better)
    off_target  : 1 − precision  — unwanted AAs encoded             (↓ better)
    f1          : harmonic mean of precision and recall              (↑ better)
    direct_obj  : 2 − recall − precision  — the value we minimise   (↓ better)
    """
    overlap    = orig_set & enc_aas
    recall     = len(overlap) / len(orig_set)  if orig_set else 1.0
    precision  = len(overlap) / len(enc_aas)   if enc_aas  else 0.0
    loss       = 1.0 - recall
    off_target = 1.0 - precision
    denom      = precision + recall
    f1         = 2.0 * precision * recall / denom if denom > 0.0 else 0.0
    direct_obj = 2.0 - recall - precision
    return dict(
        recall=recall, precision=precision, loss=loss,
        off_target=off_target, f1=f1, direct_obj=direct_obj,
    )


# ═══════════════════════════════════════════════════════════════
# Codon search
# ═══════════════════════════════════════════════════════════════

def best_codon_direct(orig_set):
    """
    Exhaustively search all IUPAC codons and return the one that minimises
    2 − recall − precision directly.  No weights, always globally optimal.

    Returns (codon_str, enc_aas_frozenset, metrics_dict).
    On ties, the codon with the fewest encoded AAs is preferred (tightest
    library that is equally good on the objective).
    """
    best_str  = None
    best_aas  = frozenset()
    best_obj  = float("inf")
    best_size = float("inf")

    for codon, enc_aas in CODON_AA_MAP.items():
        m   = position_metrics(orig_set, enc_aas)
        obj = m["direct_obj"]
        sz  = len(enc_aas)
        # Primary: minimise objective; secondary: minimise library size on tie
        if obj < best_obj or (obj == best_obj and sz < best_size):
            best_obj  = obj
            best_str  = codon
            best_aas  = enc_aas
            best_size = sz

    return best_str, best_aas, position_metrics(orig_set, best_aas)


def best_codon_weighted(orig_set, w_loss, w_offtarget):
    """
    Exhaustively search all IUPAC codons and return the one that minimises
    w_loss * loss + w_offtarget * off_target.

    Returns (codon_str, enc_aas_frozenset, metrics_dict).
    Tie-breaking: smallest library size, same as direct mode.
    """
    best_str  = None
    best_aas  = frozenset()
    best_wobj = float("inf")
    best_size = float("inf")

    for codon, enc_aas in CODON_AA_MAP.items():
        m    = position_metrics(orig_set, enc_aas)
        wobj = w_loss * m["loss"] + w_offtarget * m["off_target"]
        sz   = len(enc_aas)
        if wobj < best_wobj or (wobj == best_wobj and sz < best_size):
            best_wobj = wobj
            best_str  = codon
            best_aas  = enc_aas
            best_size = sz

    return best_str, best_aas, position_metrics(orig_set, best_aas)


# ═══════════════════════════════════════════════════════════════
# Alignment pass
# ═══════════════════════════════════════════════════════════════

def run_pass(cols, pick_fn, verbose=False):
    """
    Run one full alignment pass using pick_fn to select a codon per column.

    pick_fn(orig_set) must return (codon_str, enc_aas, metrics_dict).

    Returns
    -------
    results  : list of per-position result dicts (gap-only columns omitted)
    sequence : joined IUPAC codon string  ('NNN' inserted at gap-only columns)
    """
    results  = []
    sequence = []

    for i, col in enumerate(cols):
        orig_set = set(col) - {"-"}

        if not orig_set:
            sequence.append("NNN")
            continue

        codon, enc_aas, m = pick_fn(orig_set)
        junk = enc_aas - orig_set

        results.append(dict(
            pos      = i + 1,
            orig_set = orig_set,
            enc_aas  = enc_aas,
            junk     = junk,
            codon    = codon,
            H        = shannon_entropy(col),
            **m,
        ))
        sequence.append(codon)

        if verbose:
            print(f"  Pos {i+1}")
            print(f"    Entropy (bits)  : {results[-1]['H']:.3f}")
            print(f"    Observed AAs    : {sorted(orig_set)}")
            print(f"    Encoded AAs     : {sorted(enc_aas)}")
            print(f"    Off-target AAs  : {sorted(junk)}")
            print(f"    Codon           : {codon}")
            print(f"    Recall          : {m['recall']:.3f}")
            print(f"    Precision       : {m['precision']:.3f}")
            print(f"    Off-target rate : {m['off_target']:.3f}")
            print(f"    F1              : {m['f1']:.3f}")
            print(f"    2−R−P           : {m['direct_obj']:.4f}\n")

    return results, "".join(sequence)


def summarise(results, label=""):
    """Print a summary table for a completed pass."""
    if not results:
        print("  (no valid positions)\n")
        return

    n = len(results)
    avg = lambda k: sum(r[k] for r in results) / n

    orig_space = 1
    enc_space  = 1
    junk_total = 0
    enc_total  = 0
    for r in results:
        orig_space *= max(len(r["orig_set"]), 1)
        enc_space  *= max(len(r["enc_aas"]),  1)
        junk_total += len(r["junk"])
        enc_total  += len(r["enc_aas"])

    compression  = orig_space / enc_space if enc_space > 0 else float("inf")
    global_ot    = junk_total / enc_total  if enc_total  > 0 else 0.0

    if label:
        print(f"  [{label}]")
    print(f"  Avg recall           : {avg('recall'):.4f}  (coverage of observed AAs)")
    print(f"  Avg precision        : {avg('precision'):.4f}  (fraction of encoded AAs that were observed)")
    print(f"  Avg loss             : {avg('loss'):.4f}  (missed observed AAs)")
    print(f"  Avg off-target rate  : {avg('off_target'):.4f}  (junk AAs encoded)")
    print(f"  Avg F1               : {avg('f1'):.4f}")
    print(f"  Avg 2−R−P            : {avg('direct_obj'):.4f}  (primary objective, lower = better)")
    print(f"  Compression ratio    : {compression:.4f}  (original÷encoded variant space)")
    print(f"  Global off-target    : {global_ot:.4f}  (junk AA count ÷ total encoded AA count)")
    print()


# ═══════════════════════════════════════════════════════════════
# Grid search
# ═══════════════════════════════════════════════════════════════

def _make_weight_range(grid_max, grid_step):
    """
    Return evenly-spaced floats from 0 to grid_max (inclusive).
    Uses integer arithmetic to avoid floating-point drift
    (e.g.  0.1+0.1+0.1 ≠ 0.3 in IEEE 754).
    """
    steps = round(grid_max / grid_step)
    return [round(i * grid_step, 8) for i in range(steps + 1)]


def grid_search(cols, grid_max=2.0, grid_step=0.1):
    """
    Sweep w_loss ∈ [0, grid_max] and w_offtarget ∈ [0, grid_max] in steps
    of grid_step.  For each (w_l, w_ot) pair run a full weighted pass and
    record the resulting avg(2−recall−precision).

    Returns
    -------
    grid_results : list of (w_loss, w_offtarget, avg_direct_obj)
                   in row-major order (w_loss outer loop)
    best_weights : (w_loss, w_offtarget) that achieved the lowest avg_direct_obj
    best_obj     : float — the lowest avg_direct_obj found
    """
    ws    = _make_weight_range(grid_max, grid_step)
    total = len(ws) ** 2
    done  = 0

    grid_results = []
    best_obj     = float("inf")
    best_weights = (ws[0], ws[0])   # safe initialisation

    for w_l in ws:
        for w_ot in ws:
            pick = lambda s, wl=w_l, wot=w_ot: best_codon_weighted(s, wl, wot)
            res, _ = run_pass(cols, pick, verbose=False)
            done += 1

            if not res:
                continue

            avg_dobj = sum(r["direct_obj"] for r in res) / len(res)
            grid_results.append((w_l, w_ot, avg_dobj))

            if avg_dobj < best_obj:
                best_obj     = avg_dobj
                best_weights = (w_l, w_ot)

            # Progress indicator — overwrite the same line
            pct = 100.0 * done / total
            print(f"  Grid search: {pct:5.1f}%  ({done}/{total})", end="\r")

    print()   # move past the progress line
    return grid_results, best_weights, best_obj


# ═══════════════════════════════════════════════════════════════
# ASCII heatmap
# ═══════════════════════════════════════════════════════════════

def print_ascii_heatmap(grid_results, grid_max, grid_step):
    """
    Print a text heatmap of avg(2−recall−precision) over the weight grid.

    Cells are placed by nearest-index lookup on the weight values stored in
    each grid_results entry, so the mapping is correct even when some cells
    produced no results and were skipped by grid_search.

    Column header width = 4 chars (" 0.0").
    Data cell width     = 4 chars ("  + ").
    Both derived from the same CELL_W constant so they stay aligned.

    Shading (all visible, no space character):
        + = best (lowest obj)  →  .  →  :  →  -  →  # = worst (highest)
    """
    ws = _make_weight_range(grid_max, grid_step)
    n  = len(ws)

    # ── place each result into the correct grid cell ────────────────────────
    # We look up the nearest index by value rather than assuming a contiguous
    # flat index — this is correct even when cells are skipped.
    grid = [[None] * n for _ in range(n)]
    for w_l, w_ot, obj in grid_results:
        i_l  = min(range(n), key=lambda k: abs(ws[k] - w_l))
        i_ot = min(range(n), key=lambda k: abs(ws[k] - w_ot))
        grid[i_l][i_ot] = obj

    all_vals = [v for row in grid for v in row if v is not None]
    if not all_vals:
        print("  (no data to display)\n")
        return

    lo   = min(all_vals)
    hi   = max(all_vals)
    span = hi - lo if hi > lo else 1.0

    # 5-level ASCII shading — no space so every cell is always visible
    SHADES = "+.:-#"   # index 0 = best, index 4 = worst

    def shade(val):
        if val is None:
            return "?"
        norm = (val - lo) / span                      # 0.0 = best, 1.0 = worst
        idx  = int(norm * (len(SHADES) - 1) + 0.5)   # round to nearest level
        return SHADES[max(0, min(len(SHADES) - 1, idx))]

    CELL_W  = 4   # header and data cells are both this many characters wide
    ROW_LBL = 6   # row label field width

    print("\n  Heatmap: avg(2 − recall − precision) over weight grid")
    print("  Y-axis: w_loss  |  X-axis: w_offtarget")
    print(f"  Shading:  + = best (lowest)  →  # = worst (highest)\n")

    # Header: each column label right-justified into CELL_W chars
    header = " " * (ROW_LBL + 1)
    for w in ws:
        header += f"{w:.1f}".rjust(CELL_W)
    print(header)

    # Rows: label + one cell per column, each cell = "  {ch} " = CELL_W chars
    for i_l, w_l in enumerate(ws):
        row = f"{w_l:.2f}".rjust(ROW_LBL) + " "
        for i_ot in range(n):
            row += f"  {shade(grid[i_l][i_ot])} "
        print(row)

    print(f"\n  Range: {lo:.4f} (best '+') → {hi:.4f} (worst '#')")

    # ── show the actual distinct objective values behind the shading ────────
    # When the landscape is nearly flat (e.g. only 2 distinct values) the
    # shading alone looks binary and confusing.  This always shows the real
    # numbers so you know whether the variation is meaningful or negligible.
    unique_vals = sorted(set(round(v, 6) for v in all_vals))
    print(f"\n  Distinct avg(2−R−P) values across all {len(all_vals)} cells:")
    for v in unique_vals:
        ch    = shade(v)
        count = sum(1 for row in grid for x in row
                    if x is not None and round(x, 6) == v)
        bar   = "█" * min(count, 40)
        print(f"    '{ch}'  {v:.6f}  ({count:4d} cells)  {bar}")

    if len(unique_vals) <= 3:
        print(
            "\n  ℹ  Very few distinct values — the weight landscape is essentially\n"
            "     flat once w_loss > 0.  The w_loss=0 row is a special case\n"
            "     (recall is not penalised at all, so the optimiser ignores it).\n"
            "     All other weight combinations converge to the same optimum,\n"
            "     confirming that --mode direct already finds the global solution."
        )


# ═══════════════════════════════════════════════════════════════
# Polynomial response surface  (degree 2, no numpy)
# ═══════════════════════════════════════════════════════════════

def fit_surface(grid_results):
    """
    Fit a degree-2 polynomial response surface:
        z = a + b·x + c·y + d·x² + e·x·y + f·y²
    where x = w_loss, y = w_offtarget, z = avg_direct_obj.

    Solved via ordinary least squares (normal equations + Gaussian
    elimination with partial pivoting).  No external dependencies.

    Returns
    -------
    coeffs   : (a, b, c, d, e, f)
    minimum  : (w_loss, w_offtarget) — analytical minimum of the surface;
               may lie outside the grid if the surface is not bowl-shaped
    z_min    : predicted objective at the analytical minimum
    """
    p = 6   # number of polynomial terms

    # Build design matrix rows and target values
    X_rows = []
    z_vals = []
    for x, y, z in grid_results:
        X_rows.append([1.0, x, y, x*x, x*y, y*y])
        z_vals.append(z)

    # Accumulate XᵀX and Xᵀz
    XtX = [[0.0]*p for _ in range(p)]
    Xtz = [0.0]*p
    for row, z in zip(X_rows, z_vals):
        for i in range(p):
            Xtz[i] += row[i] * z
            for j in range(p):
                XtX[i][j] += row[i] * row[j]

    # Solve via Gauss-Jordan elimination with partial pivoting
    aug = [XtX[i][:] + [Xtz[i]] for i in range(p)]
    for col in range(p):
        # Find pivot
        max_row = max(range(col, p), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-12:
            continue
        # Normalise pivot row
        aug[col] = [v / pivot for v in aug[col]]
        # Eliminate column in all other rows
        for r in range(p):
            if r == col:
                continue
            factor = aug[r][col]
            aug[r] = [aug[r][c] - factor * aug[col][c] for c in range(p + 1)]

    coeffs = [aug[i][p] for i in range(p)]
    a, b, c, d, e, f = coeffs

    # Analytical minimum: solve ∂z/∂x = 0, ∂z/∂y = 0
    #   2d·x + e·y = −b
    #   e·x + 2f·y = −c
    det = 4.0*d*f - e*e
    if abs(det) < 1e-12:
        # Degenerate surface (flat or saddle) — return grid minimum
        best = min(grid_results, key=lambda t: t[2])
        return coeffs, (best[0], best[1]), best[2]

    x_min = (-2.0*f*b + e*c) / det
    y_min = (-2.0*d*c + e*b) / det
    z_min = a + b*x_min + c*y_min + d*x_min**2 + e*x_min*y_min + f*y_min**2

    return coeffs, (x_min, y_min), z_min


def surface_r_squared(grid_results, coeffs):
    """Coefficient of determination (R²) for the fitted surface."""
    a, b, c, d, e, f = coeffs
    z_mean = sum(t[2] for t in grid_results) / len(grid_results)
    ss_res = ss_tot = 0.0
    for x, y, z in grid_results:
        z_pred = a + b*x + c*y + d*x*x + e*x*y + f*y*y
        ss_res += (z - z_pred)**2
        ss_tot += (z - z_mean)**2
    return 1.0 - ss_res/ss_tot if ss_tot > 1e-12 else 1.0


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="Degenerate IUPAC codon design from a protein MSA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("fasta",
        help="Aligned FASTA file of protein sequences")
    p.add_argument("--mode", choices=["direct","grid","surface"],
        default="direct",
        help="Optimisation mode (default: direct)")
    p.add_argument("--grid-step", type=float, default=0.1, metavar="F",
        help="Step size for weight grid (default: 0.1)")
    p.add_argument("--grid-max",  type=float, default=2.0, metavar="F",
        help="Upper bound for weight grid (default: 2.0)")
    p.add_argument("--w-loss", type=float, default=None, metavar="F",
        help="Manual weight for recall loss (must be >= 0). "
             "Supply together with --w-offtarget to run a single weighted "
             "pass instead of a grid/surface sweep.")
    p.add_argument("--w-offtarget", type=float, default=None, metavar="F",
        help="Manual weight for off-target rate (must be >= 0). "
             "Supply together with --w-loss.")
    p.add_argument("--verbose", action="store_true",
        help="Print per-position detail during the final pass")
    args = p.parse_args()

    # Validate: both weights must be given together, and must be non-negative
    if (args.w_loss is None) != (args.w_offtarget is None):
        p.error("--w-loss and --w-offtarget must be supplied together.")
    if args.w_loss is not None and args.w_loss < 0:
        p.error("--w-loss must be >= 0.")
    if args.w_offtarget is not None and args.w_offtarget < 0:
        p.error("--w-offtarget must be >= 0.")

    return args


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

SEP = "=" * 58

def main():
    args = parse_args()
    seqs = read_fasta(args.fasta)
    cols = columns(seqs)

    # ── Direct ─────────────────────────────────────────────────
    # ── Manual weights (overrides --mode) ──────────────────────
    if args.w_loss is not None:
        wl, wot = args.w_loss, args.w_offtarget
        print(SEP)
        print(f"MODE: manual weights  —  w_loss={wl}  w_offtarget={wot}")
        print(SEP + "\n")

        pick = lambda s: best_codon_weighted(s, wl, wot)
        results, sequence = run_pass(cols, pick, verbose=args.verbose)

        # Also run direct to show how far off the global optimum we are
        res_direct, _ = run_pass(cols, best_codon_direct, verbose=False)
        direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
        manual_obj = sum(r["direct_obj"] for r in results)    / len(results)
        gap        = manual_obj - direct_obj

        print(SEP)
        print("RESULTS")
        print(SEP + "\n")
        print(f"Degenerate sequence:\n  {sequence}\n")
        summarise(results, label=f"w_loss={wl}, w_offtarget={wot}")
        print(f"  Direct optimum avg(2-R-P) : {direct_obj:.6f}")
        print(f"  This run   avg(2-R-P)     : {manual_obj:.6f}")
        verdict = "identical to direct" if gap < 1e-4 else f"{gap:+.6f} above direct optimum"
        print(f"  Gap                        : {verdict}\n")
        return

    if args.mode == "direct":
        print(SEP)
        print("MODE: direct  —  minimise 2 − recall − precision")
        print(SEP + "\n")

        results, sequence = run_pass(cols, best_codon_direct, verbose=args.verbose)

        print(SEP)
        print("RESULTS")
        print(SEP + "\n")
        print(f"Degenerate sequence:\n  {sequence}\n")
        summarise(results)

    # ── Grid ───────────────────────────────────────────────────
    elif args.mode == "grid":
        print(SEP)
        print("MODE: grid search")
        print(SEP + "\n")

        # Direct optimum — gold standard
        res_direct, seq_direct = run_pass(cols, best_codon_direct, verbose=False)
        direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
        print(f"Direct optimum  avg(2−R−P) = {direct_obj:.6f}\n")

        n_combos = round((args.grid_max / args.grid_step + 1) ** 2)
        print(f"Running {n_combos} weight combinations...\n")

        grid_results, best_w, best_grid_obj = grid_search(
            cols, args.grid_max, args.grid_step
        )

        print_ascii_heatmap(grid_results, args.grid_max, args.grid_step)

        print(f"\nBest grid weights     : w_loss={best_w[0]:.2f},  "
              f"w_offtarget={best_w[1]:.2f}")
        print(f"Best grid avg(2−R−P)  : {best_grid_obj:.6f}")
        print(f"Direct avg(2−R−P)     : {direct_obj:.6f}")
        gap = best_grid_obj - direct_obj
        verdict = "≈ identical to direct" if gap < 1e-4 else "grid is suboptimal"
        print(f"Gap                   : {gap:.6f}  ({verdict})\n")

        # Final pass at best weights
        pick = lambda s: best_codon_weighted(s, best_w[0], best_w[1])
        res_grid, seq_grid = run_pass(cols, pick, verbose=args.verbose)

        print(SEP)
        print("RESULTS")
        print(SEP + "\n")
        print(f"Degenerate sequence (best grid weights):\n  {seq_grid}\n")
        summarise(res_grid, label=f"w_loss={best_w[0]:.2f}, w_offtarget={best_w[1]:.2f}")

        print(f"Degenerate sequence (direct):\n  {seq_direct}\n")
        summarise(res_direct, label="direct")

    # ── Surface ────────────────────────────────────────────────
    elif args.mode == "surface":
        print(SEP)
        print("MODE: polynomial response surface (degree 2)")
        print(SEP + "\n")

        # Direct optimum — gold standard
        res_direct, seq_direct = run_pass(cols, best_codon_direct, verbose=False)
        direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
        print(f"Direct optimum  avg(2−R−P) = {direct_obj:.6f}\n")

        n_combos = round((args.grid_max / args.grid_step + 1) ** 2)
        print(f"Running {n_combos} weight combinations...\n")

        grid_results, best_w_grid, best_grid_obj = grid_search(
            cols, args.grid_max, args.grid_step
        )

        print_ascii_heatmap(grid_results, args.grid_max, args.grid_step)

        # Fit and report surface
        coeffs, (x_min, y_min), z_min = fit_surface(grid_results)
        r2 = surface_r_squared(grid_results, coeffs)
        a, b, c, d, e, f = coeffs

        print(f"\nPolynomial surface fit  (R² = {r2:.4f})")
        print(f"  z = {a:+.4f}"
              f" {b:+.4f}·wl"
              f" {c:+.4f}·wot"
              f" {d:+.4f}·wl²"
              f" {e:+.4f}·wl·wot"
              f" {f:+.4f}·wot²")

        in_range = (0.0 <= x_min <= args.grid_max and
                    0.0 <= y_min <= args.grid_max)
        print(f"\nAnalytical minimum of surface:")
        print(f"  w_loss      = {x_min:.4f}")
        print(f"  w_offtarget = {y_min:.4f}")
        print(f"  predicted z = {z_min:.6f}")
        if not in_range:
            print(f"  ⚠  Outside [0, {args.grid_max}] — optimum is at the boundary.")

        print(f"\nBest grid weights     : w_loss={best_w_grid[0]:.2f},  "
              f"w_offtarget={best_w_grid[1]:.2f}")
        print(f"Best grid avg(2−R−P)  : {best_grid_obj:.6f}")
        print(f"Direct avg(2−R−P)     : {direct_obj:.6f}")

        # Run with surface-suggested weights (clipped to valid range)
        x_use = max(0.0, min(args.grid_max, x_min))
        y_use = max(0.0, min(args.grid_max, y_min))
        pick  = lambda s: best_codon_weighted(s, x_use, y_use)
        res_surf, seq_surf = run_pass(cols, pick, verbose=args.verbose)
        surf_actual = sum(r["direct_obj"] for r in res_surf) / len(res_surf)

        gap = surf_actual - direct_obj
        verdict = "≈ identical to direct" if gap < 1e-4 else "surface is suboptimal"
        print(f"\nActual avg(2−R−P) at surface weights : {surf_actual:.6f}")
        print(f"Gap vs direct                         : {gap:.6f}  ({verdict})\n")

        print(SEP)
        print("RESULTS")
        print(SEP + "\n")
        print(f"Degenerate sequence (surface weights):\n  {seq_surf}\n")
        summarise(res_surf, label=f"w_loss={x_use:.4f}, w_offtarget={y_use:.4f}")

        print(f"Degenerate sequence (direct):\n  {seq_direct}\n")
        summarise(res_direct, label="direct")

        print(SEP)
        print("INTERPRETING THE SURFACE FIT")
        print(SEP)
        print("""
  R² ≈ 1.0  The objective landscape is a smooth quadratic bowl.
            The analytical minimum is a reliable guide to optimal
            weights.  surface and direct sequences should match.

  R² < 0.8  The landscape is irregular or nearly flat — weight
            tuning makes little difference.  Use direct mode.

  Analytical minimum outside [0, grid-max]
            The optimal weight is 0 for one or both penalties,
            meaning the optimiser should ignore that term entirely.
            Direct mode already achieves this implicitly.

  Direct mode always produces the globally optimal sequence for
  the 2−recall−precision objective, regardless of surface quality.
  Grid/surface modes are diagnostic tools to understand the weight
  landscape and verify that the weighted proxy objective agrees
  with the direct objective.
""")


if __name__ == "__main__":
    main()