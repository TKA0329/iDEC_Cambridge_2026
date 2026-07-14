"""
reverse_translation_V4_Binned.py
=================================
Extends V3 with Hamming-distance-based clustering (k-medoids) before
reverse translation.  Instead of producing one degenerate sequence from
the full sequence list, V4 splits the list into k bins, runs reverse
translation on each bin, and reports the total encoded sequence space
across all bins — which is typically smaller (tighter) than the single
degenerate sequence produced from the full list.

All three V3 optimisation modes (direct, grid, surface) are available
and apply per-bin.  Manual weight override (--w-loss / --w-offtarget)
is also supported.  The --bins-sweep feature always uses direct mode
(grid/surface per-bin in a sweep would be prohibitively slow).

CLI options
-----------
Binning (new in V4):
    --bins K          Cluster input into K bins and produce K degenerate
                      sequences.  Default: 1 (reproduces V3 behaviour).
    --bins-sweep      Sweep K = 1..5 in direct mode and print a
                      comparison table.  Overrides --bins and --mode.
    --sample-size N   Randomly sample N sequences before clustering.
                      Default: 1000.  Set 0 to use the full list
                      (warning: very slow for large inputs).
    --seed S          Random seed for sampling and clustering (default: 42).

Optimisation (inherited from V3 — apply per-bin):
    --mode            direct | grid | surface   (default: direct)
    --grid-step F     Step size for grid/surface sweep (default: 0.1)
    --grid-max  F     Upper bound of weight sweep (default: 2.0)
    --w-loss F        Manual recall-loss weight (supply with --w-offtarget)
    --w-offtarget F   Manual off-target weight  (supply with --w-loss)
    --verbose         Print per-position detail for each bin

Usage examples
--------------
    # Single run, 3 bins, direct mode (default)
    python reverse_translation_V4_Binned.py sequences.fasta --bins 3

    # 3 bins with grid optimisation per bin
    python reverse_translation_V4_Binned.py sequences.fasta --bins 3 --mode grid

    # 3 bins with manual weights per bin
    python reverse_translation_V4_Binned.py sequences.fasta --bins 3 \\
        --w-loss 0.5 --w-offtarget 1.2

    # Sweep k=1..5 (always direct mode)
    python reverse_translation_V4_Binned.py sequences.fasta --bins-sweep

    # Larger sample, fixed seed
    python reverse_translation_V4_Binned.py sequences.fasta --bins-sweep \\
        --sample-size 2000 --seed 7
"""

from collections import Counter, defaultdict
from itertools import product as iproduct
import argparse
import math
import random
import time
import sys


# ═══════════════════════════════════════════════════════════════
# Progress bar
# ═══════════════════════════════════════════════════════════════

class ProgressBar:
    """
    Single-line terminal progress bar with ETA.

    Usage
    -----
        pb = ProgressBar(total=n, label="Clustering")
        for i, item in enumerate(items):
            pb.update(i + 1)
            ...do work...
        pb.finish()
    """
    BAR_WIDTH = 30

    def __init__(self, total, label=""):
        self.total     = max(total, 1)
        self.label     = label
        self.start     = time.perf_counter()
        self._last_len = 0   # chars written on current line — for clean overwrite
        self.update(0)

    def update(self, done):
        elapsed  = time.perf_counter() - self.start
        fraction = done / self.total
        filled   = int(fraction * self.BAR_WIDTH)
        bar      = "█" * filled + "░" * (self.BAR_WIDTH - filled)
        pct      = fraction * 100

        if done > 0 and elapsed > 0:
            rate     = done / elapsed
            remaining = (self.total - done) / rate
            eta_str  = _fmt_seconds(remaining)
        else:
            eta_str  = "--:--"

        elapsed_str = _fmt_seconds(elapsed)
        line = (
            f"  {self.label}  [{bar}] "
            f"{pct:5.1f}%  {done}/{self.total}  "
            f"elapsed {elapsed_str}  ETA {eta_str}"
        )
        # Overwrite previous line
        sys.stdout.write("\r" + line)
        sys.stdout.flush()
        self._last_len = len(line)

    def finish(self):
        """Print completed bar on its own line."""
        self.update(self.total)
        sys.stdout.write("\n")
        sys.stdout.flush()


def _fmt_seconds(s):
    """Format a duration in seconds as  Xh Ym Zs  or  Mm Zs  or  Zs."""
    s = int(s)
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    if h:
        return f"{h}h {m:02d}m {sec:02d}s"
    if m:
        return f"{m}m {sec:02d}s"
    return f"{sec}s"

# ═══════════════════════════════════════════════════════════════
# Tee — mirror stdout to a log file
# ═══════════════════════════════════════════════════════════════

class Tee:
    """
    Wraps sys.stdout so every write goes to both the terminal and a file.
    Progress-bar \r lines are written verbatim to the file, then a newline
    is appended so the log is readable as plain text.

    Usage
    -----
        tee = Tee("run.log")
        sys.stdout = tee
        ...
        tee.close()          # restores sys.stdout and closes the file
    """
    def __init__(self, path):
        self._terminal = sys.stdout
        self._file     = open(path, "w", encoding="utf-8")
        sys.stdout     = self

    def write(self, text):
        self._terminal.write(text)
        self._terminal.flush()
        if "\r" in text:
            # Convert carriage-return lines to newline-terminated for the log
            self._file.write(text.replace("\r", "\n"))
        else:
            self._file.write(text)
        self._file.flush()

    def flush(self):
        self._terminal.flush()
        self._file.flush()

    def close(self):
        sys.stdout = self._terminal
        self._file.close()


# ═══════════════════════════════════════════════════════════════
# Constants  (identical to V3)
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

IUPAC = {
    "A":{"A"},"C":{"C"},"G":{"G"},"T":{"T"},
    "R":{"A","G"},"Y":{"C","T"},"S":{"G","C"},
    "W":{"A","T"},"K":{"G","T"},"M":{"A","C"},
    "B":{"C","G","T"},"D":{"A","G","T"},
    "H":{"A","C","T"},"V":{"A","C","G"},"N":{"A","C","G","T"},
}


def _expand(codon):
    return ("".join(t) for t in iproduct(*[IUPAC[b] for b in codon]))


def _build_codon_aa_map():
    mapping = {}
    for c in ("".join(t) for t in iproduct(IUPAC.keys(), repeat=3)):
        aas = frozenset(
            CODON_TABLE[x]
            for x in _expand(c)
            if x in CODON_TABLE and CODON_TABLE[x] != "*"
        )
        if aas:
            mapping[c] = aas
    return mapping


CODON_AA_MAP = _build_codon_aa_map()


# ═══════════════════════════════════════════════════════════════
# I/O
# ═══════════════════════════════════════════════════════════════

def read_fasta(path):
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
    L = len(seqs[0])
    return [[s[i] for s in seqs if s[i] != "-"] for i in range(L)]


# ═══════════════════════════════════════════════════════════════
# Metrics  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def shannon_entropy(col):
    counts = Counter(col)
    total  = len(col)
    return -sum(
        (v / total) * math.log2(v / total)
        for v in counts.values() if v > 0
    )


def position_metrics(orig_set, enc_aas):
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
# Codon search  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def best_codon_direct(orig_set):
    best_str  = None
    best_aas  = frozenset()
    best_obj  = float("inf")
    best_size = float("inf")

    for codon, enc_aas in CODON_AA_MAP.items():
        m   = position_metrics(orig_set, enc_aas)
        obj = m["direct_obj"]
        sz  = len(enc_aas)
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
    Tie-breaking: smallest library size.
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
# Grid search  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def _make_weight_range(grid_max, grid_step):
    steps = round(grid_max / grid_step)
    return [round(i * grid_step, 8) for i in range(steps + 1)]


def grid_search(cols, grid_max=2.0, grid_step=0.1):
    ws    = _make_weight_range(grid_max, grid_step)
    total = len(ws) ** 2
    done  = 0

    grid_results = []
    best_obj     = float("inf")
    best_weights = (ws[0], ws[0])

    for w_l in ws:
        for w_ot in ws:
            pick = lambda s, wl=w_l, wot=w_ot: best_codon_weighted(s, wl, wot)
            res, _ = run_pass(cols, pick, verbose=False, label=f"grid {done+1:>4}/{total}")
            done += 1

            if not res:
                continue

            avg_dobj = sum(r["direct_obj"] for r in res) / len(res)
            grid_results.append((w_l, w_ot, avg_dobj))

            if avg_dobj < best_obj:
                best_obj     = avg_dobj
                best_weights = (w_l, w_ot)

    return grid_results, best_weights, best_obj


# ═══════════════════════════════════════════════════════════════
# ASCII heatmap  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def print_ascii_heatmap(grid_results, grid_max, grid_step):
    ws = _make_weight_range(grid_max, grid_step)
    n  = len(ws)

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
    SHADES = "+.:-#"

    def shade(val):
        if val is None:
            return "?"
        norm = (val - lo) / span
        idx  = int(norm * (len(SHADES) - 1) + 0.5)
        return SHADES[max(0, min(len(SHADES) - 1, idx))]

    CELL_W  = 4
    ROW_LBL = 6

    print("\n  Heatmap: avg(2 − recall − precision) over weight grid")
    print("  Y-axis: w_loss  |  X-axis: w_offtarget")
    print(f"  Shading:  + = best (lowest)  →  # = worst (highest)\n")

    header = " " * (ROW_LBL + 1)
    for w in ws:
        header += f"{w:.1f}".rjust(CELL_W)
    print(header)

    for i_l, w_l in enumerate(ws):
        row = f"{w_l:.2f}".rjust(ROW_LBL) + " "
        for i_ot in range(n):
            row += f"  {shade(grid[i_l][i_ot])} "
        print(row)

    print(f"\n  Range: {lo:.4f} (best '+') → {hi:.4f} (worst '#')")

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
            "     flat.  --mode direct already finds the global solution."
        )


# ═══════════════════════════════════════════════════════════════
# Polynomial response surface  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def fit_surface(grid_results):
    p = 6
    X_rows = []
    z_vals = []
    for x, y, z in grid_results:
        X_rows.append([1.0, x, y, x*x, x*y, y*y])
        z_vals.append(z)

    XtX = [[0.0]*p for _ in range(p)]
    Xtz = [0.0]*p
    for row, z in zip(X_rows, z_vals):
        for i in range(p):
            Xtz[i] += row[i] * z
            for j in range(p):
                XtX[i][j] += row[i] * row[j]

    aug = [XtX[i][:] + [Xtz[i]] for i in range(p)]
    for col in range(p):
        max_row = max(range(col, p), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-12:
            continue
        aug[col] = [v / pivot for v in aug[col]]
        for r in range(p):
            if r == col:
                continue
            factor = aug[r][col]
            aug[r] = [aug[r][c] - factor * aug[col][c] for c in range(p + 1)]

    coeffs = [aug[i][p] for i in range(p)]
    a, b, c, d, e, f = coeffs

    det = 4.0*d*f - e*e
    if abs(det) < 1e-12:
        best = min(grid_results, key=lambda t: t[2])
        return coeffs, (best[0], best[1]), best[2]

    x_min = (-2.0*f*b + e*c) / det
    y_min = (-2.0*d*c + e*b) / det
    z_min = a + b*x_min + c*y_min + d*x_min**2 + e*x_min*y_min + f*y_min**2
    return coeffs, (x_min, y_min), z_min


def surface_r_squared(grid_results, coeffs):
    a, b, c, d, e, f = coeffs
    z_mean = sum(t[2] for t in grid_results) / len(grid_results)
    ss_res = ss_tot = 0.0
    for x, y, z in grid_results:
        z_pred = a + b*x + c*y + d*x*x + e*x*y + f*y*y
        ss_res += (z - z_pred)**2
        ss_tot += (z - z_mean)**2
    return 1.0 - ss_res/ss_tot if ss_tot > 1e-12 else 1.0

def run_pass(cols, pick_fn, verbose=False, label="rev-translate", verbose_file=None):
    """
    verbose=True  + verbose_file=None  -> print per-position detail to stdout (k=1 behaviour).
    verbose=True  + verbose_file=<fh>  -> write per-position detail to that file handle instead.
    verbose=False                      -> no per-position detail regardless of verbose_file.
    """
    results  = []
    sequence = []
    n_cols   = len(cols)
    pb       = ProgressBar(n_cols, label=f"  {label:<14}")

    for i, col in enumerate(cols):
        orig_set = set(col) - {"-"}

        if not orig_set:
            sequence.append("NNN")
            pb.update(i + 1)
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
        pb.update(i + 1)

    pb.finish()

    if verbose:
        out = verbose_file if verbose_file is not None else sys.stdout
        for r in results:
            out.write(f"  Pos {r['pos']}\n")
            out.write(f"    Observed AAs    : {sorted(r['orig_set'])}\n")
            out.write(f"    Encoded AAs     : {sorted(r['enc_aas'])}\n")
            out.write(f"    Codon           : {r['codon']}\n")
            out.write(f"    Recall          : {r['recall']:.3f}\n")
            out.write(f"    Precision       : {r['precision']:.3f}\n")
            out.write(f"    2-R-P           : {r['direct_obj']:.4f}\n\n")

    return results, "".join(sequence)


def encoded_space(results):
    """Product of encoded AA set sizes across all positions."""
    space = 1
    for r in results:
        space *= max(len(r["enc_aas"]), 1)
    return space


def seq_coverage(seqs, results):
    """
    Fraction of input sequences fully contained within the degenerate
    codon space.  A sequence is 'covered' if, at every position that has
    a valid result, its amino acid appears in enc_aas for that position.

    Returns (n_covered, total, fraction).
    """
    # Build a position-indexed lookup: pos (1-based) -> enc_aas
    enc_at = {r["pos"]: r["enc_aas"] for r in results}

    n_covered = 0
    for seq in seqs:
        covered = True
        for r in results:
            aa = seq[r["pos"] - 1] if r["pos"] - 1 < len(seq) else "-"
            if aa == "-":
                continue
            if aa not in enc_at[r["pos"]]:
                covered = False
                break
        if covered:
            n_covered += 1

    total = len(seqs)
    fraction = n_covered / total if total > 0 else 0.0
    return n_covered, total, fraction


def summarise(results, label="", show_sequence=True, sequence=None, baseline_space=None, coverage_seqs=None):
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

    compression = orig_space / enc_space if enc_space > 0 else float("inf")
    global_ot   = junk_total / enc_total  if enc_total  > 0 else 0.0

    if label:
        print(f"  [{label}]")
    if show_sequence and sequence:
        print(f"  Degenerate sequence:\n    {sequence}\n")
    print(f"  Avg recall           : {avg('recall'):.4f}")
    print(f"  Avg precision        : {avg('precision'):.4f}")
    print(f"  Avg loss             : {avg('loss'):.4f}")
    print(f"  Avg off-target rate  : {avg('off_target'):.4f}")
    print(f"  Avg F1               : {avg('f1'):.4f}")
    print(f"  Avg 2−R−P            : {avg('direct_obj'):.4f}  (lower = better)")
    print(f"  Encoded space        : {enc_space}")
    if baseline_space is not None and baseline_space > 0:
        pct_of_baseline = enc_space / baseline_space * 100
        reduction_vs_baseline = 100.0 - pct_of_baseline
        print(f"  vs k=1 baseline      : {pct_of_baseline:.2f}% of baseline  "
              f"({reduction_vs_baseline:+.2f}% reduction)")
    print(f"  Compression ratio    : {compression:.4f}  (original÷encoded)")
    print(f"  Global off-target    : {global_ot:.4f}")
    if coverage_seqs is not None:
        n_cov, n_tot, frac = seq_coverage(coverage_seqs, results)
        print(f"  Seqs covered         : {n_cov}/{n_tot}  ({frac*100:.2f}%  of overall sequence set fully within degenerate space)")
    print()
    return enc_space


# ═══════════════════════════════════════════════════════════════
# Hamming clustering  (k-medoids / PAM)
# ═══════════════════════════════════════════════════════════════

def hamming_distance(s1, s2):
    """Number of positions at which two aligned sequences differ."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def cluster_sequences(sequences, k, max_iter=100, seed=42):
    """
    K-medoids clustering using Hamming distance.

    Parameters
    ----------
    sequences : list of str  — aligned protein sequences (equal length)
    k         : int          — number of clusters
    max_iter  : int          — maximum PAM iterations
    seed      : int          — random seed for medoid initialisation

    Returns
    -------
    List of k lists, each containing the sequences assigned to that cluster.
    Empty clusters are dropped silently.
    """
    if k >= len(sequences):
        # Degenerate case: more bins than sequences
        return [[s] for s in sequences]

    rng = random.Random(seed)
    medoids = rng.sample(sequences, k)

    n = len(sequences)

    for iteration in range(max_iter):
        # Assignment step
        clusters = defaultdict(list)
        pb_assign = ProgressBar(n, label=f"iter {iteration+1:>3}  assign ")
        for idx, seq in enumerate(sequences):
            closest = min(medoids, key=lambda m: hamming_distance(seq, m))
            clusters[closest].append(seq)
            pb_assign.update(idx + 1)
        pb_assign.finish()

        # Update step — work unit = n_cluster^2 per cluster
        total_update_work = sum(
            len(clusters.get(m, [m])) ** 2 for m in medoids
        )
        pb_update = ProgressBar(total_update_work, label=f"iter {iteration+1:>3}  update")
        work_done = 0
        new_medoids = []
        for medoid in medoids:
            cluster_seqs = clusters.get(medoid, [medoid])
            csz = len(cluster_seqs)
            best = min(
                cluster_seqs,
                key=lambda s1: sum(hamming_distance(s1, s2) for s2 in cluster_seqs)
            )
            new_medoids.append(best)
            work_done += csz * csz
            pb_update.update(work_done)
        pb_update.finish()

        if set(new_medoids) == set(medoids):
            print(f"  Clustering converged after {iteration + 1} iteration(s).")
            break
        medoids = new_medoids
    else:
        print(f"  Clustering reached max_iter={max_iter} without convergence.")

    return [clusters[m] for m in medoids if clusters.get(m)]


# ═══════════════════════════════════════════════════════════════
# Binned reverse translation
# ═══════════════════════════════════════════════════════════════

def run_binned(seqs, k, sample_size, seed, pick_fn, verbose=False,
               overall_pb=None, overall_done=None, overall_total=None,
               verbose_fh=None):
    """
    Cluster seqs into k bins, run reverse translation (using pick_fn) on
    each bin, and return (bin_results, total_encoded_space).

    pick_fn:    codon-selection function — best_codon_direct or a
                functools.partial wrapping best_codon_weighted.
    overall_pb / overall_done: optional shared ProgressBar for sweep mode.
    verbose_fh: open file handle for per-position verbose output.
                When provided and verbose=True and k>1, each bin's detail
                is written there (with a section header) instead of stdout.
                The caller is responsible for opening and closing the file.

    Returns (bin_results, total_encoded_space, working) where each element
    of bin_results is (results, sequence, bin_size, bin_seqs), and working
    is the full sampled-or-full sequence list used for this k (i.e. the
    union across all bins) — the right set to test overall coverage against.
    """
    n_orig = len(seqs)

    # Sampling
    if sample_size > 0 and sample_size < n_orig:
        rng = random.Random(seed)
        working = rng.sample(seqs, sample_size)
        print(f"  Sampled {sample_size} / {n_orig} sequences for clustering.")
    else:
        working = seqs
        print(f"  Using all {n_orig} sequences (no sampling).")

    if k == 1:
        print(f"  k=1: skipping clustering, running on full sample.")
        cols = columns(working)
        # k=1: verbose stays on stdout as in V3
        # k=1: verbose stays on stdout as in V3 (no file needed)
        results, sequence = run_pass(cols, pick_fn, verbose=verbose, label="bin 1 / 1   ")
        total_space = encoded_space(results)
        if overall_pb is not None:
            overall_pb.update(overall_done[0] + 1)
            overall_done[0] += 1
        return [(results, sequence, len(working), working)], total_space, working

    # Cluster
    print(f"  Clustering into k={k} bins...")
    bins = cluster_sequences(working, k=k, seed=seed)
    actual_k = len(bins)
    if actual_k < k:
        print(f"  Warning: only {actual_k} non-empty bins produced (requested {k}).")

    # Reverse translate each bin
    bin_results = []
    total_space = 0
    for i, bin_seqs in enumerate(bins, 1):
        print(f"  Bin {i}/{actual_k}  ({len(bin_seqs)} sequences):")
        cols = columns(bin_seqs)

        # Verbose for binned runs: write to shared file handle with section header
        vfile_to_use = None
        if verbose and verbose_fh is not None:
            vfile_to_use = verbose_fh
            verbose_fh.write(f"\n{'=' * 60}\n")
            verbose_fh.write(f"BIN {i}/{actual_k}  ({len(bin_seqs)} sequences)\n")
            verbose_fh.write(f"{'=' * 60}\n\n")

        results, sequence = run_pass(
            cols, pick_fn, verbose=verbose,
            label=f"bin {i}/{actual_k}      ",
            verbose_file=vfile_to_use,
        )

        space = encoded_space(results)
        total_space += space
        bin_results.append((results, sequence, len(bin_seqs), bin_seqs))
        if overall_pb is not None:
            # Flush a newline so the overall bar \r doesn't stomp the previous line
            sys.stdout.write("\n")
            sys.stdout.flush()
            overall_pb.update(overall_done[0] + 1)
            overall_done[0] += 1

    return bin_results, total_space, working


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

SEP = "=" * 60

def parse_args():
    p = argparse.ArgumentParser(
        description="Binned degenerate IUPAC codon design (V4).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("fasta",
        help="Aligned FASTA file of protein sequences")

    # ── Binning (V4) ───────────────────────────────────────────
    p.add_argument("--bins", type=int, default=1, metavar="K",
        help="Number of bins/clusters (default: 1 = V3 behaviour)")
    p.add_argument("--bins-sweep", action="store_true",
        help="Sweep k=1..N in direct mode and print comparison table "
             "(overrides --bins and --mode)")
    p.add_argument("--sweep-max", type=int, default=5, metavar="N",
        help="Maximum k to sweep up to when using --bins-sweep (default: 5)")
    p.add_argument("--sample-size", type=int, default=1000, metavar="N",
        help="Sequences to sample before clustering (default: 1000; 0 = no sampling)")
    p.add_argument("--seed", type=int, default=42, metavar="S",
        help="Random seed for sampling and clustering (default: 42)")

    # ── Optimisation mode (V3) — apply per-bin ─────────────────
    p.add_argument("--mode", choices=["direct", "grid", "surface"],
        default="direct",
        help="Optimisation mode per bin (default: direct). "
             "Ignored when --bins-sweep is set.")
    p.add_argument("--grid-step", type=float, default=0.1, metavar="F",
        help="Step size for weight grid (default: 0.1)")
    p.add_argument("--grid-max",  type=float, default=2.0, metavar="F",
        help="Upper bound for weight grid (default: 2.0)")
    p.add_argument("--w-loss", type=float, default=None, metavar="F",
        help="Manual recall-loss weight (supply together with --w-offtarget)")
    p.add_argument("--w-offtarget", type=float, default=None, metavar="F",
        help="Manual off-target weight (supply together with --w-loss)")

    p.add_argument("--verbose", action="store_true",
        help="Print per-position detail for each bin")
    p.add_argument("--csv", type=str, default=None, metavar="FILE",
        help="Write per-bin summary table to a CSV file "
             "(k, bin, bin_size, encoded_space, pct_of_baseline, reduction_pct, "
             "seqs_covered_pct, degenerate_sequence)")
    p.add_argument("--log", type=str, default=None, metavar="FILE",
        help="Mirror all terminal output to FILE for future reference")
    args = p.parse_args()

    # Validate manual weights and sweep-max
    if args.sweep_max < 1:
        p.error("--sweep-max must be >= 1.")
    if (args.w_loss is None) != (args.w_offtarget is None):
        p.error("--w-loss and --w-offtarget must be supplied together.")
    if args.w_loss is not None and args.w_loss < 0:
        p.error("--w-loss must be >= 0.")
    if args.w_offtarget is not None and args.w_offtarget < 0:
        p.error("--w-offtarget must be >= 0.")

    return args


def main():
    args = parse_args()

    # ── Start log tee before any output ────────────────────────
    tee = None
    if args.log:
        tee = Tee(args.log)
        # Write a header so the log is self-contained
        import datetime
        print(f"# reverse_translation_V4_Binned.py  log")
        print(f"# Date : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"# Args : {' '.join(sys.argv[1:])}")
        print(f"# {'=' * 56}")

    seqs = read_fasta(args.fasta)
    print(f"\nLoaded {len(seqs)} sequences, length {len(seqs[0])} AA each.")

    # ── Resolve pick function ───────────────────────────────────
    # Manual weights override --mode entirely (same as V3).
    # In sweep mode, always use direct regardless of --mode.
    if args.w_loss is not None:
        wl, wot = args.w_loss, args.w_offtarget
        pick_fn  = lambda s: best_codon_weighted(s, wl, wot)
        mode_str = f"manual weights  w_loss={wl}  w_offtarget={wot}"
    else:
        pick_fn  = best_codon_direct
        mode_str = "direct"   # overridden below for grid/surface in single-run path

    # ── Open single verbose output file if requested ──────────
    verbose_fh = None
    if args.verbose:
        import os
        fasta_base = os.path.splitext(os.path.basename(args.fasta))[0]
        verbose_path = f"verbose_{fasta_base}.txt"
        verbose_fh = open(verbose_path, "w")
        verbose_fh.write(f"Verbose per-position detail\n")
        verbose_fh.write(f"Input : {args.fasta}\n")
        verbose_fh.write(f"Mode  : {'bins-sweep' if args.bins_sweep else args.mode}\n")
        verbose_fh.write(f"Bins  : {'sweep k=1..5' if args.bins_sweep else args.bins}\n")
        verbose_fh.write("=" * 60 + "\n")
        print(f"  Verbose detail will be written to: {verbose_path}")

    if args.bins_sweep:
        # ── Sweep k=1..5 (always direct) ───────────────────────
        print(SEP)
        sweep_max = args.sweep_max
        print(f"MODE: bins sweep  —  k = 1 .. {sweep_max}  (direct per bin)")
        print(f"Sample size per run : {args.sample_size if args.sample_size > 0 else 'all'}")
        print(SEP + "\n")

        sweep_max      = args.sweep_max
        sweep_results  = []
        baseline_space = None
        total_passes   = sweep_max * (sweep_max + 1) // 2   # 1+2+...+sweep_max
        overall_pb     = ProgressBar(total_passes, label="overall progress")
        overall_done   = [0]
        csv_rows       = []   # accumulated for optional CSV output

        for k in range(1, sweep_max + 1):
            print(SEP)
            print(f"k = {k}")
            print(SEP)
            bin_results, total_space, working = run_binned(
                seqs, k=k,
                sample_size=args.sample_size,
                seed=args.seed,
                pick_fn=best_codon_direct,
                verbose=args.verbose,
                overall_pb=overall_pb,
                overall_done=overall_done,
                overall_total=total_passes,
                verbose_fh=verbose_fh if args.verbose and k > 1 else None,
            )
            if k == 1:
                baseline_space = total_space
            reduction = (1.0 - total_space / baseline_space) * 100 if baseline_space else 0.0
            sweep_results.append((k, total_space, len(bin_results), reduction))

            sys.stdout.write("\n")
            sys.stdout.flush()
            for i, (results, sequence, bin_size, bin_seqs) in enumerate(bin_results, 1):
                bin_space = encoded_space(results)
                n_cov, n_tot, frac = seq_coverage(working, results)
                if baseline_space and baseline_space > 0:
                    pct = bin_space / baseline_space * 100
                    reduction_bin = 100.0 - pct
                    space_str = f"{bin_space}  ({pct:.2f}% of k=1 baseline  /  {reduction_bin:.2f}% reduction)"
                else:
                    pct          = 100.0
                    reduction_bin = 0.0
                    space_str    = str(bin_space)
                print(f"  Bin {i} ({bin_size} sequences):")
                print(f"    Degenerate sequence: {sequence}")
                print(f"    Encoded space      : {space_str}")
                print(f"    Seqs covered       : {n_cov}/{n_tot}  ({frac*100:.2f}%)")
                csv_rows.append({
                    "k":                    k,
                    "bin":                  i,
                    "bin_size":             bin_size,
                    "encoded_space":        bin_space,
                    "pct_of_baseline":      round(pct, 4),
                    "reduction_pct":        round(reduction_bin, 4),
                    "seqs_covered_pct":     round(frac * 100, 4),
                    "degenerate_sequence":  sequence,
                })
            print()

        overall_pb.finish()

        print(SEP)
        print("COMPARISON TABLE")
        print(SEP)
        print(f"  {'k':>4}  {'Total enc. space':>20}  {'Oligos':>6}  {'Reduction vs k=1':>18}")
        print(f"  {'-'*4}  {'-'*20}  {'-'*6}  {'-'*18}")
        for k, space, n_oligos, reduction in sweep_results:
            marker = "  ← baseline" if k == 1 else ""
            print(f"  {k:>4}  {space:>20}  {n_oligos:>6}  {reduction:>17.2f}%{marker}")
        print()
        print("  Interpretation:")
        print("  - Total encoded space = sum of variant spaces across all bin oligos.")
        print("  - Lower space = tighter library = less off-target coverage.")
        print("  - More bins = more oligos to order; weigh efficiency vs synthesis cost.")
        print()

    else:
        # ── Single k run ────────────────────────────────────────
        k = args.bins

        # For grid/surface, we need cols — resolve after clustering.
        # The pick_fn for grid/surface is determined per-bin inside run_binned_with_mode.
        if args.w_loss is not None:
            mode_label = f"manual weights  w_loss={args.w_loss}  w_offtarget={args.w_offtarget}"
        else:
            mode_label = args.mode

        print(SEP)
        print(f"MODE: {mode_label}  —  k = {k}")
        print(f"Sample size : {args.sample_size if args.sample_size > 0 else 'all'}")
        print(SEP + "\n")

        # k=1 baseline always uses direct for fair comparison
        print("Running k=1 baseline (direct)...")
        baseline_results, total_space_baseline, baseline_working = run_binned(
            seqs, k=1,
            sample_size=args.sample_size,
            seed=args.seed,
            pick_fn=best_codon_direct,
            verbose=False,
        )
        baseline_seq, baseline_bin_seqs = baseline_results[0][1], baseline_results[0][3]
        print()

        if args.mode in ("grid", "surface") and args.w_loss is None:
            # Grid/surface: run per-bin with grid search
            # We need to sample + cluster first, then run grid on each bin's cols
            print(f"Running k={k} with {args.mode} mode per bin...")
            n_orig  = len(seqs)
            rng     = random.Random(args.seed)
            working = rng.sample(seqs, args.sample_size) if (args.sample_size > 0 and args.sample_size < n_orig) else seqs
            print(f"  {'Sampled' if args.sample_size > 0 and args.sample_size < n_orig else 'Using all'} "
                  f"{len(working)} sequences.")

            if k == 1:
                bins = [working]
            else:
                print(f"  Clustering into k={k} bins...")
                bins = cluster_sequences(working, k=k, seed=args.seed)

            bin_results = []
            total_space_k = 0

            for i, bin_seqs in enumerate(bins, 1):
                actual_k = len(bins)
                print(f"\n{SEP}")
                print(f"  Bin {i}/{actual_k}  ({len(bin_seqs)} sequences) — {args.mode} mode")
                print(SEP)
                cols = columns(bin_seqs)

                # Direct baseline for this bin
                res_direct, seq_direct = run_pass(cols, best_codon_direct, verbose=False,
                                                   label=f"bin{i} direct  ")
                direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
                print(f"  Direct optimum  avg(2−R−P) = {direct_obj:.6f}\n")

                n_combos = round((args.grid_max / args.grid_step + 1) ** 2)
                print(f"  Running {n_combos} weight combinations...")
                grid_res, best_w, best_grid_obj = grid_search(cols, args.grid_max, args.grid_step)
                print_ascii_heatmap(grid_res, args.grid_max, args.grid_step)

                # Verbose for grid/surface binned runs: write to shared verbose_fh
                vfile_to_use = None
                if args.verbose and verbose_fh is not None and k > 1:
                    vfile_to_use = verbose_fh
                    verbose_fh.write(f"\n{'=' * 60}\n")
                    verbose_fh.write(f"k={k}  {args.mode.upper()} MODE  BIN {i}/{actual_k}"
                                     f"  ({len(bin_seqs)} sequences)\n")
                    verbose_fh.write(f"{'=' * 60}\n\n")

                if args.mode == "grid":
                    gap = best_grid_obj - direct_obj
                    verdict = "≈ identical to direct" if gap < 1e-4 else "grid is suboptimal"
                    print(f"\n  Best grid weights : w_loss={best_w[0]:.2f}  w_offtarget={best_w[1]:.2f}")
                    print(f"  Best grid avg(2−R−P) : {best_grid_obj:.6f}  ({verdict})\n")
                    pick = lambda s, wl=best_w[0], wot=best_w[1]: best_codon_weighted(s, wl, wot)
                    results, sequence = run_pass(cols, pick, verbose=args.verbose,
                                                  label=f"bin{i} grid    ",
                                                  verbose_file=vfile_to_use)
                    summarise(results, label=f"bin {i} (grid best weights)", sequence=sequence,
                              baseline_space=total_space_baseline, coverage_seqs=working)

                else:  # surface
                    coeffs, (x_min, y_min), z_min = fit_surface(grid_res)
                    r2 = surface_r_squared(grid_res, coeffs)
                    a, b, c, d, e, f = coeffs
                    in_range = (0.0 <= x_min <= args.grid_max and 0.0 <= y_min <= args.grid_max)
                    print(f"\n  Surface fit  R² = {r2:.4f}")
                    print(f"  Analytical minimum:  w_loss={x_min:.4f}  w_offtarget={y_min:.4f}")
                    if not in_range:
                        print(f"  ⚠  Outside [0, {args.grid_max}] — optimum at boundary.")
                    x_use = max(0.0, min(args.grid_max, x_min))
                    y_use = max(0.0, min(args.grid_max, y_min))
                    pick  = lambda s, wl=x_use, wot=y_use: best_codon_weighted(s, wl, wot)
                    results, sequence = run_pass(cols, pick, verbose=args.verbose,
                                                  label=f"bin{i} surface ",
                                                  verbose_file=vfile_to_use)
                    summarise(results, label=f"bin {i} (surface weights)", sequence=sequence,
                              baseline_space=total_space_baseline, coverage_seqs=working)

                space = encoded_space(results)
                total_space_k += space
                bin_results.append((results, sequence, len(bin_seqs), bin_seqs))

        else:
            # Direct or manual weights — use run_binned as normal
            print(f"Running k={k}...")
            bin_results, total_space_k, working_k = run_binned(
                seqs, k=k,
                sample_size=args.sample_size,
                seed=args.seed,
                pick_fn=pick_fn,
                verbose=args.verbose,
                verbose_fh=verbose_fh if args.verbose and k > 1 else None,
            )

            print(SEP)
            print(f"RESULTS  (k={k})")
            print(SEP + "\n")

            print(f"  Baseline (k=1, direct):")
            print(f"    Degenerate sequence : {baseline_seq}")
            print(f"    Encoded space       : {total_space_baseline}\n")

            for i, (results, sequence, bin_size, bin_seqs) in enumerate(bin_results, 1):
                print(f"  Bin {i} ({bin_size} sequences):")
                summarise(results, label=f"bin {i}", sequence=sequence,
                          baseline_space=total_space_baseline, coverage_seqs=working_k)

            # Also show direct comparison if manual weights were used
            if args.w_loss is not None:
                res_direct, _ = run_pass(columns(seqs[:min(len(seqs), args.sample_size or len(seqs))]),
                                          best_codon_direct, verbose=False, label="direct check  ")
                direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
                manual_obj = sum(r["direct_obj"] for r in bin_results[0][0]) / len(bin_results[0][0])
                gap = manual_obj - direct_obj
                verdict = "identical to direct" if gap < 1e-4 else f"{gap:+.6f} above direct optimum"
                print(f"  Direct optimum avg(2-R-P) : {direct_obj:.6f}")
                print(f"  Manual weights avg(2-R-P) : {manual_obj:.6f}  ({verdict})\n")

        reduction = (1.0 - total_space_k / total_space_baseline) * 100
        print(SEP)
        print("SUMMARY")
        print(SEP)
        print(f"  k=1 encoded space (direct baseline) : {total_space_baseline}")
        print(f"  k={k} total encoded space            : {total_space_k}")
        print(f"  Space reduction                      : {reduction:.2f}%")
        print(f"  Oligos to order                      : {len(bin_results)}  (vs 1 for k=1)")
        print()


    # ── Write CSV if requested ─────────────────────────────────
    if args.csv and args.bins_sweep:
        import csv
        fieldnames = ["k", "bin", "bin_size", "encoded_space",
                      "pct_of_baseline", "reduction_pct", "seqs_covered_pct",
                      "degenerate_sequence"]
        with open(args.csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"  CSV written to: {args.csv}")
    elif args.csv and not args.bins_sweep:
        print("  Note: --csv is only populated during --bins-sweep runs.")

    # ── Close verbose file ─────────────────────────────────────
    if verbose_fh is not None:
        verbose_fh.close()
        print(f"\n  Verbose detail written to: {verbose_path}")

    # ── Close log tee ───────────────────────────────────────────
    if tee is not None:
        tee.close()
        # Print to real stdout now that tee is closed
        print(f"  Terminal output written to: {args.log}")


if __name__ == "__main__":
    main()