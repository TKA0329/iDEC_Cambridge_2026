"""
reverse_translation_V5_StopSafe.py
===================================
Extends V4 with stop-codon awareness.

The problem: an IUPAC degenerate codon like "NNN" expands to 64 explicit
triplets. Some of those triplets (TAA, TAG, TGA) are stop codons. V4's
CODON_AA_MAP silently *dropped* stop-codon expansions when computing
which amino acids a degenerate codon encodes — which made the objective
function (recall/precision) look correct, but ignored the fact that a
real oligo pool synthesized from that degenerate codon will, some
fraction of the time, actually synthesize a premature stop at that
position. That's a synthesis risk V4 never surfaced.

V5 fixes this by building two separate codon search spaces:

    CODON_AA_MAP_SAFE  — only degenerate codons whose FULL expansion is
                          100% free of stop codons. No triplet the
                          codon can collapse to is TAA/TAG/TGA. This is
                          the new DEFAULT search space, so a synthesized
                          oligo can never truncate early.

    CODON_AA_MAP_ALL   — the old V4 behaviour: every degenerate codon
                          that encodes at least one non-stop amino acid,
                          even if some of its expansions are stop codons.
                          Available via --allow-stop-codons for users who
                          want the (larger, more compressive) V4 search
                          space and are willing to accept some stop-codon
                          risk in the synthesized pool.

New stats
---------
Default (stop-safe) mode prints a "cost of stop-safety" comparison:
how much worse (if at all) the stop-free codon choice is, in terms of
avg loss / off-target rate, versus the old unrestricted V4 search per
position.

--allow-stop-codons mode prints an estimated percentage of the encoded
sequence space that contains AT LEAST ONE spontaneously formed stop
codon — i.e. the fraction of synthesized oligos from that degenerate
sequence that would be expected to truncate early, assuming uniform
independent base draws at each degenerate position.

Everything else (binning, k-medoids clustering, direct/grid/surface
weight optimisation, bins-sweep, CSV export, verbose logging) is
inherited unchanged from V4.

Bug fix (V5): fully-gapped alignment columns used to be hardcoded to
the codon "NNN", which carries stop-codon risk (3/64 expansions) and
silently ignored whichever codon_map was active. They now go through
pick_fn like any other column, targeting the full 20-aa alphabet, so
they respect --allow-stop-codons and stay stop-free by default.

CLI options
-----------
Stop-codon handling (new in V5):
    --allow-stop-codons   Search the full V4 IUPAC codon space (may
                           spontaneously encode stop codons at some
                           fraction of synthesized oligos). Default is
                           OFF — only stop-free codons are considered.

Binning (from V4):
    --bins K          Cluster input into K bins and produce K degenerate
                      sequences.  Default: 1 (reproduces V3 behaviour).
    --bins-sweep      Sweep K = 1..N in direct mode and print a
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
    # Single run, 3 bins, direct mode, stop-safe (default)
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins 3

    # Same, but allow the old V4 unrestricted codon space
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins 3 \\
        --allow-stop-codons

    # 3 bins with grid optimisation per bin
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins 3 --mode grid

    # 3 bins with manual weights per bin
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins 3 \\
        --w-loss 0.5 --w-offtarget 1.2

    # Sweep k=1..5 (always direct mode)
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins-sweep

    # Larger sample, fixed seed
    python reverse_translation_V5_StopSafe.py sequences.fasta --bins-sweep \\
        --sample-size 2000 --seed 7
"""

from collections import Counter, defaultdict
from itertools import product as iproduct
import argparse
import math
import random
import time
import sys

# NumPy is optional. When present, clustering uses vectorized Hamming
# distances (large speedup on --sample-size 0 / big inputs). When absent,
# everything falls back to the original pure-Python implementation.
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    np = None
    NUMPY_AVAILABLE = False


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


def _build_codon_maps():
    """
    Build the codon search spaces used by best_codon_direct /
    best_codon_weighted.

    CODON_AA_MAP_ALL   : every degenerate codon that encodes at least one
                          non-stop amino acid (V4 behaviour). The amino
                          acid set for each codon excludes stop-codon
                          expansions, but the codon itself may still
                          expand to TAA/TAG/TGA some fraction of the
                          time — CODON_STOP_FRAC records that fraction.

    CODON_AA_MAP_SAFE  : the subset of CODON_AA_MAP_ALL whose full
                          expansion is 100% free of stop codons (i.e.
                          CODON_STOP_FRAC[codon] == 0.0). This is the
                          new V5 default search space, guaranteeing a
                          synthesized oligo can never terminate early
                          at that position.

    CODON_STOP_FRAC    : codon -> fraction of the codon's full explicit
                          expansion (all triplets it can collapse to,
                          counted with multiplicity) that translate to
                          a stop codon. Only defined for codons that
                          appear in CODON_AA_MAP_ALL (i.e. codons whose
                          expansion is not 100% stop codons).
    """
    aa_map_all = {}
    stop_frac  = {}
    for c in ("".join(t) for t in iproduct(IUPAC.keys(), repeat=3)):
        expansions = list(_expand(c))
        total = len(expansions)
        stops = sum(1 for x in expansions if CODON_TABLE.get(x) == "*")
        aas = frozenset(
            CODON_TABLE[x] for x in expansions
            if x in CODON_TABLE and CODON_TABLE[x] != "*"
        )
        if aas:
            aa_map_all[c] = aas
            stop_frac[c]  = stops / total
    aa_map_safe = {c: aas for c, aas in aa_map_all.items() if stop_frac[c] == 0.0}
    return aa_map_all, aa_map_safe, stop_frac


CODON_AA_MAP_ALL, CODON_AA_MAP_SAFE, CODON_STOP_FRAC = _build_codon_maps()

# Per-codon penalty used by best_codon_stop_aware: -log(1 - stop_frac).
# Summing this across positions equals -log(sequence-level stop-free
# probability), i.e. it's the exact additive decomposition of the true
# multiplicative risk computed by stop_risk_estimate() -- not a linear
# approximation of it. Safe for every codon in CODON_AA_MAP_ALL since
# none of them have stop_frac == 1.0 (each has at least one non-stop AA
# by construction), so this never hits log(0).
CODON_STOP_LOG_PENALTY = {c: -math.log(1.0 - p) for c, p in CODON_STOP_FRAC.items()}

# Backwards-compatible default: the safe map, since it's now the default
# search space almost everywhere in this file.
CODON_AA_MAP = CODON_AA_MAP_SAFE

# The full 20-amino-acid alphabet, used as a "no information" target set
# for fully-gapped columns (see run_pass) so the wildcard codon chosen
# for those positions is still selected from the active codon_map rather
# than being hardcoded to "NNN" (which has nonzero stop-codon risk).
ALL_AMINO_ACIDS = frozenset(aa for aa in CODON_TABLE.values() if aa != "*")


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

def best_codon_direct(orig_set, codon_map=CODON_AA_MAP_SAFE):
    """
    codon_map defaults to CODON_AA_MAP_SAFE (stop-free codons only).
    Pass codon_map=CODON_AA_MAP_ALL to reproduce V4's unrestricted search.
    """
    best_str  = None
    best_aas  = frozenset()
    best_obj  = float("inf")
    best_size = float("inf")

    for codon, enc_aas in codon_map.items():
        m   = position_metrics(orig_set, enc_aas)
        obj = m["direct_obj"]
        sz  = len(enc_aas)
        if obj < best_obj or (obj == best_obj and sz < best_size):
            best_obj  = obj
            best_str  = codon
            best_aas  = enc_aas
            best_size = sz

    return best_str, best_aas, position_metrics(orig_set, best_aas)


def best_codon_weighted(orig_set, w_loss, w_offtarget, codon_map=CODON_AA_MAP_SAFE):
    """
    Exhaustively search codon_map and return the codon that minimises
    w_loss * loss + w_offtarget * off_target.
    Tie-breaking: smallest library size.

    codon_map defaults to CODON_AA_MAP_SAFE (stop-free codons only).
    Pass codon_map=CODON_AA_MAP_ALL to reproduce V4's unrestricted search.
    """
    best_str  = None
    best_aas  = frozenset()
    best_wobj = float("inf")
    best_size = float("inf")

    for codon, enc_aas in codon_map.items():
        m    = position_metrics(orig_set, enc_aas)
        wobj = w_loss * m["loss"] + w_offtarget * m["off_target"]
        sz   = len(enc_aas)
        if wobj < best_wobj or (wobj == best_wobj and sz < best_size):
            best_wobj = wobj
            best_str  = codon
            best_aas  = enc_aas
            best_size = sz

    return best_str, best_aas, position_metrics(orig_set, best_aas)


def best_codon_stop_aware(orig_set, w_loss, w_offtarget, w_stop, codon_map=CODON_AA_MAP_ALL):
    """
    Like best_codon_weighted, but adds a third weighted term penalizing
    stop-codon risk:

        w_loss * loss + w_offtarget * off_target + w_stop * stop_penalty

    codon_map defaults to CODON_AA_MAP_ALL (the unrestricted space), not
    CODON_AA_MAP_SAFE -- every codon in the safe map already has zero
    stop risk, so w_stop would have nothing to act on there. Pass
    CODON_AA_MAP_SAFE explicitly if you want w_stop to just be inert.

    stop_penalty is -log(1 - CODON_STOP_FRAC[codon]), not the raw
    stop_frac. This isn't cosmetic: the quantity you actually care about
    -- the probability a synthesized oligo has at least one spontaneous
    stop, as computed by stop_risk_estimate() -- is
    1 - product(1 - p_i) across positions, which is multiplicative.
    Summing raw stop_frac per position only approximates that risk when
    every p_i is small, and understates it once any position carries
    real risk. Summing -log(1 - p_i) instead makes the position-by-
    position sum exactly equal to -log(true sequence-level survival
    probability), so minimizing it here is an exact decomposition of
    minimizing the true sequence-level risk, not a first-order
    approximation of it.

    w_stop=0 recovers best_codon_weighted's ranking over codon_map.
    Increasing w_stop trades recall/precision for lower stop risk; a
    large w_stop effectively forces stop-free codons wherever one
    exists in codon_map for that position.

    Tie-breaking: smallest library size (same as best_codon_weighted).
    """
    best_str  = None
    best_aas  = frozenset()
    best_wobj = float("inf")
    best_size = float("inf")

    for codon, enc_aas in codon_map.items():
        m        = position_metrics(orig_set, enc_aas)
        stop_pen = CODON_STOP_LOG_PENALTY.get(codon, 0.0)
        wobj     = w_loss * m["loss"] + w_offtarget * m["off_target"] + w_stop * stop_pen
        sz       = len(enc_aas)
        if wobj < best_wobj or (wobj == best_wobj and sz < best_size):
            best_wobj = wobj
            best_str  = codon
            best_aas  = enc_aas
            best_size = sz

    return best_str, best_aas, position_metrics(orig_set, best_aas)


# ═══════════════════════════════════════════════════════════════
# Stop-codon risk stats  (new in V5)
# ═══════════════════════════════════════════════════════════════

def stop_risk_estimate(results):
    """
    Estimate the probability that a randomly synthesized oligo drawn
    from this degenerate sequence contains AT LEAST ONE spontaneously
    formed in-frame stop codon.

    Assumes each IUPAC base at each degenerate position is drawn
    independently and uniformly at random from its allowed base set
    (the standard assumption for a mixed-base oligo pool). Under that
    assumption, the probability the whole oligo is stop-free is the
    product across positions of (1 - per-position stop fraction), so
    P(>=1 stop) = 1 - that product.

    Returns (p_at_least_one_stop, per_position) where per_position is
    a list of (pos, codon, stop_fraction_at_that_position).
    """
    survival = 1.0
    per_position = []
    for r in results:
        codon  = r["codon"]
        p_stop = CODON_STOP_FRAC.get(codon, 0.0)
        per_position.append((r["pos"], codon, p_stop))
        survival *= (1.0 - p_stop)
    return 1.0 - survival, per_position


def compare_codon_maps(cols):
    """
    Run direct-mode codon selection twice on the same columns — once
    restricted to CODON_AA_MAP_SAFE (stop-free, the V5 default) and
    once over CODON_AA_MAP_ALL (unrestricted, V4 behaviour) — and
    return average loss / off-target / objective for both, plus the
    deltas. This quantifies the "cost of stop-safety": how much worse
    (if at all) restricting the search to stop-free codons is,
    compared to the old unrestricted V4 search.
    """
    res_safe, _ = run_pass(
        cols, lambda s: best_codon_direct(s, CODON_AA_MAP_SAFE),
        verbose=False, label="stop-safe check",
    )
    res_all, _ = run_pass(
        cols, lambda s: best_codon_direct(s, CODON_AA_MAP_ALL),
        verbose=False, label="all-codon check ",
    )
    avg = lambda res, key: (sum(r[key] for r in res) / len(res)) if res else 0.0
    return dict(
        safe_loss=avg(res_safe, "loss"),
        safe_offtarget=avg(res_safe, "off_target"),
        safe_obj=avg(res_safe, "direct_obj"),
        all_loss=avg(res_all, "loss"),
        all_offtarget=avg(res_all, "off_target"),
        all_obj=avg(res_all, "direct_obj"),
    )


def stop_codon_summary_line(cols, results, allow_stop_codons):
    """One-line version of print_stop_codon_stats, for dense per-bin output
    (e.g. inside --bins-sweep, where a full breakdown per bin would be
    too noisy)."""
    if not allow_stop_codons:
        cmp = compare_codon_maps(cols)
        d_loss = cmp["all_loss"] - cmp["safe_loss"]
        d_ot   = cmp["all_offtarget"] - cmp["safe_offtarget"]
        return f"stop-safety cost vs unrestricted: Δloss {d_loss:+.4f}  Δoff-target {d_ot:+.4f}"
    else:
        p_stop, _ = stop_risk_estimate(results)
        return f"est. sequence space with >=1 spontaneous stop: {p_stop * 100:.2f}%"


def print_stop_codon_stats(cols, results, allow_stop_codons, label=""):
    """
    Print the V5 stop-codon stat appropriate to the active mode.

    Stop-safe mode (default): print the safe-vs-all comparison from
    compare_codon_maps — the cost, in avg loss/off-target, of
    restricting the search to codons that can never spontaneously
    form a stop codon.

    --allow-stop-codons mode: print the estimated percentage of the
    encoded sequence space (i.e. fraction of synthesized oligos) that
    contains at least one spontaneous stop codon, from
    stop_risk_estimate.
    """
    tag = f" [{label}]" if label else ""
    if not allow_stop_codons:
        cmp = compare_codon_maps(cols)
        d_loss = cmp["all_loss"] - cmp["safe_loss"]
        d_ot   = cmp["all_offtarget"] - cmp["safe_offtarget"]
        d_obj  = cmp["all_obj"] - cmp["safe_obj"]
        print(f"  Stop-safety cost{tag}  (stop-free search vs unrestricted V4 search):")
        print(f"    Avg loss        : safe {cmp['safe_loss']:.4f}  vs  all {cmp['all_loss']:.4f}  "
              f"(Δ {d_loss:+.4f})")
        print(f"    Avg off-target  : safe {cmp['safe_offtarget']:.4f}  vs  all {cmp['all_offtarget']:.4f}  "
              f"(Δ {d_ot:+.4f})")
        print(f"    Avg 2−R−P       : safe {cmp['safe_obj']:.4f}  vs  all {cmp['all_obj']:.4f}  "
              f"(Δ {d_obj:+.4f})")
        if abs(d_loss) < 1e-6 and abs(d_ot) < 1e-6:
            print(f"    -> No cost: the stop-free codon set achieves the same optimum here.")
        else:
            print(f"    -> Restricting to stop-free codons costs the above, in exchange for "
                  f"guaranteeing no premature stop.")
        print()
        return cmp
    else:
        p_stop, per_position = stop_risk_estimate(results)
        risky = [p for p in per_position if p[2] > 0.0]
        print(f"  Spontaneous stop-codon risk{tag}  (--allow-stop-codons mode):")
        print(f"    Est. % of sequence space with >=1 stop codon : {p_stop * 100:.4f}%")
        print(f"    Positions using a codon with nonzero stop risk : {len(risky)}/{len(per_position)}")
        if risky:
            worst = sorted(risky, key=lambda p: -p[2])[:5]
            print(f"    Highest-risk positions (pos, codon, per-position stop fraction):")
            for pos, codon, frac in worst:
                print(f"      pos {pos:>4}  codon {codon}  {frac*100:.2f}% of that position's draws")
        print()
        return p_stop


# ═══════════════════════════════════════════════════════════════
# Grid search  (identical to V3)
# ═══════════════════════════════════════════════════════════════

def _make_weight_range(grid_max, grid_step):
    steps = round(grid_max / grid_step)
    return [round(i * grid_step, 8) for i in range(steps + 1)]


def _grid_search_worker(task):
    """
    Module-level (picklable) worker for one grid cell. Defined at top level
    rather than as a closure/lambda so it works under both 'fork' and
    'spawn' multiprocessing start methods. Runs run_pass with
    show_progress=False: concurrent worker processes writing \r-updated
    progress bars to stdout would just interleave and corrupt each
    other's output, so progress is tracked by a single bar in the main
    process instead (one tick per completed cell).
    """
    cols, w_l, w_ot, codon_map = task
    pick = lambda s, wl=w_l, wot=w_ot: best_codon_weighted(s, wl, wot, codon_map)
    res, _ = run_pass(cols, pick, verbose=False, show_progress=False)
    if not res:
        return (w_l, w_ot, None)
    avg_dobj = sum(r["direct_obj"] for r in res) / len(res)
    return (w_l, w_ot, avg_dobj)


def grid_search(cols, grid_max=2.0, grid_step=0.1, codon_map=CODON_AA_MAP_SAFE, n_jobs=1):
    """
    Every (w_loss, w_offtarget) grid cell is fully independent — a full
    run_pass over all columns for one fixed weight pair — so this is
    embarrassingly parallel. With n_jobs > 1, cells are farmed out across
    a multiprocessing Pool; n_jobs=1 (default) keeps the original
    single-process behaviour.
    """
    ws    = _make_weight_range(grid_max, grid_step)
    total = len(ws) ** 2
    tasks = [(cols, w_l, w_ot, codon_map) for w_l in ws for w_ot in ws]

    grid_results = []
    best_obj     = float("inf")
    best_weights = (ws[0], ws[0])

    pb   = ProgressBar(total, label="  grid search   ")
    done = 0

    if n_jobs and n_jobs > 1:
        import multiprocessing as mp
        with mp.Pool(processes=n_jobs) as pool:
            for w_l, w_ot, avg_dobj in pool.imap_unordered(_grid_search_worker, tasks):
                done += 1
                pb.update(done)
                if avg_dobj is None:
                    continue
                grid_results.append((w_l, w_ot, avg_dobj))
                if avg_dobj < best_obj:
                    best_obj     = avg_dobj
                    best_weights = (w_l, w_ot)
    else:
        for task in tasks:
            w_l, w_ot, avg_dobj = _grid_search_worker(task)
            done += 1
            pb.update(done)
            if avg_dobj is None:
                continue
            grid_results.append((w_l, w_ot, avg_dobj))
            if avg_dobj < best_obj:
                best_obj     = avg_dobj
                best_weights = (w_l, w_ot)

    pb.finish()
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

def run_pass(cols, pick_fn, verbose=False, label="rev-translate", verbose_file=None,
             show_progress=True):
    """
    verbose=True  + verbose_file=None  -> print per-position detail to stdout (k=1 behaviour).
    verbose=True  + verbose_file=<fh>  -> write per-position detail to that file handle instead.
    verbose=False                      -> no per-position detail regardless of verbose_file.
    show_progress=False                -> skip the progress bar entirely. Use this for calls
                                           made from parallel workers (e.g. grid_search with
                                           --n-jobs > 1), since concurrent processes writing
                                           \r-updated bars to stdout would just interleave and
                                           corrupt each other's output.
    """
    results  = []
    sequence = []
    n_cols   = len(cols)
    pb       = ProgressBar(n_cols, label=f"  {label:<14}") if show_progress else None

    for i, col in enumerate(cols):
        orig_set = set(col) - {"-"}

        if not orig_set:
            # Fully-gapped column: no observed amino acids to encode.
            # Previously hardcoded to "NNN", which has nonzero stop-codon
            # risk (3/64 of its expansions are stop codons) and silently
            # bypassed whichever codon_map pick_fn is using. Instead, ask
            # pick_fn for its broadest-coverage codon against the full
            # 20-aa alphabet — this stays inside the active codon_map, so
            # stop-safe mode still guarantees zero stop-codon risk here,
            # and --allow-stop-codons still gets V4's old NNN-equivalent
            # behaviour if that's genuinely the best choice in that map.
            wildcard_codon, _, _ = pick_fn(ALL_AMINO_ACIDS)
            sequence.append(wildcard_codon)
            if pb is not None:
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
        if pb is not None:
            pb.update(i + 1)

    if pb is not None:
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


def _encode_matrix(seqs):
    """Encode a list of equal-length strings as an (n, L) uint8 NumPy array."""
    return np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(len(seqs), -1)


def _assign_clusters_numpy(seq_matrix, sequences, medoids):
    """
    Vectorized assignment step: for each sequence, find its closest medoid.
    Replaces the O(n*k*L) pure-Python zip/generator inner loop with k
    vectorized array comparisons (one per medoid), each O(n*L) in
    compiled code rather than the interpreter.

    Ties are broken the same way as the original `min(medoids, key=...)`
    call: the first medoid (in `medoids` list order) achieving the
    minimum distance wins, since np.argmin also returns the first
    occurrence of the minimum along an axis.
    """
    medoid_matrix = _encode_matrix(medoids)
    dist_matrix = np.empty((seq_matrix.shape[0], len(medoids)), dtype=np.int32)
    for j in range(len(medoids)):
        dist_matrix[:, j] = (seq_matrix != medoid_matrix[j]).sum(axis=1)
    closest_idx = dist_matrix.argmin(axis=1)

    clusters = defaultdict(list)
    for i, seq in enumerate(sequences):
        clusters[medoids[closest_idx[i]]].append(seq)
    return clusters


def _best_medoid_numpy(cluster_seqs):
    """Vectorized O(m^2) pairwise Hamming search for the point in
    cluster_seqs minimizing total distance to every other point in it."""
    arr = _encode_matrix(cluster_seqs)
    diffs = arr[:, None, :] != arr[None, :, :]
    totals = diffs.sum(axis=(1, 2))
    return cluster_seqs[int(totals.argmin())]


def _best_medoid_python(cluster_seqs):
    """Pure-Python fallback for the same O(m^2) pairwise search."""
    return min(
        cluster_seqs,
        key=lambda s1: sum(hamming_distance(s1, s2) for s2 in cluster_seqs)
    )


def cluster_sequences(sequences, k, max_iter=100, seed=42, medoid_cap=2000,
                       show_progress=True):
    """
    K-medoids clustering using Hamming distance (PAM).

    Parameters
    ----------
    sequences  : list of str — aligned protein sequences (equal length)
    k          : int         — number of clusters
    max_iter   : int         — maximum PAM iterations
    seed       : int         — random seed for medoid init + subsampling
    medoid_cap : int         — CLARA-style cap on cluster size used when
                 recomputing a medoid each iteration. The medoid-update
                 step is normally O(cluster_size^2), which is what makes
                 --sample-size 0 slow on large inputs: doubling a
                 cluster's size roughly quadruples this step. When a
                 cluster exceeds medoid_cap, a random subsample of that
                 size is used to choose the new medoid instead of the
                 full cluster, bounding the update step to
                 O(medoid_cap^2) regardless of how large the input gets.
                 The chosen medoid is still validated against — and the
                 resulting cluster still contains — every sequence
                 actually assigned to it; only the *search* for the best
                 medoid is subsampled. Set 0 to disable (use the full
                 cluster every time, i.e. the original exact behaviour —
                 slow for large bins).
    show_progress : bool — print a short per-iteration status line.

    Uses NumPy-vectorized Hamming distances when NumPy is available
    (~50-100x faster in practice), falling back to a pure-Python
    implementation otherwise. In both cases medoid_cap bounds the
    quadratic update step, which is the change that actually lets
    --sample-size 0 scale rather than just running the same O(n^2)
    algorithm faster.

    Returns
    -------
    List of k lists, each containing the sequences assigned to that
    cluster. Empty clusters are dropped silently.
    """
    if k >= len(sequences):
        # Degenerate case: more bins than sequences
        return [[s] for s in sequences]

    rng = random.Random(seed)
    medoids = rng.sample(sequences, k)
    n = len(sequences)
    use_numpy = NUMPY_AVAILABLE
    cap = medoid_cap if (medoid_cap and medoid_cap > 0) else None

    seq_matrix = _encode_matrix(sequences) if use_numpy else None
    clusters = None

    for iteration in range(max_iter):
        # ── Assignment step ──────────────────────────────────────
        if use_numpy:
            clusters = _assign_clusters_numpy(seq_matrix, sequences, medoids)
        else:
            clusters = defaultdict(list)
            pb_assign = ProgressBar(n, label=f"iter {iteration+1:>3}  assign ") if show_progress else None
            for idx, seq in enumerate(sequences):
                closest = min(medoids, key=lambda m: hamming_distance(seq, m))
                clusters[closest].append(seq)
                if pb_assign is not None:
                    pb_assign.update(idx + 1)
            if pb_assign is not None:
                pb_assign.finish()

        if show_progress:
            tag = "vectorized" if use_numpy else "python"
            print(f"  iter {iteration+1:>3}  assign  ({tag}, {n} sequences)")

        # ── Update step (CLARA-capped) ───────────────────────────
        # Track new_clusters keyed by the *new* medoid, built in lockstep
        # with new_medoids, so the pairing below is always internally
        # consistent — this also fixes a bug in the previous version:
        # when clustering hit max_iter without converging, `medoids` got
        # reassigned to new_medoids but `clusters` stayed keyed by the
        # *old* medoids, so the final `clusters[m] for m in medoids`
        # lookup could come back empty. Now clusters and medoids are
        # always updated together, so that mismatch can't happen.
        new_medoids = []
        new_clusters = {}
        for medoid in medoids:
            cluster_seqs = clusters.get(medoid, [medoid])
            if cap is not None and len(cluster_seqs) > cap:
                candidate_seqs = rng.sample(cluster_seqs, cap)
            else:
                candidate_seqs = cluster_seqs

            if use_numpy:
                best = _best_medoid_numpy(candidate_seqs)
            else:
                best = _best_medoid_python(candidate_seqs)

            new_medoids.append(best)
            new_clusters[best] = cluster_seqs   # full membership, keyed by new medoid

        if show_progress:
            capped = any(len(clusters.get(m, [m])) > (cap or float("inf")) for m in medoids)
            print(f"  iter {iteration+1:>3}  update  "
                  f"({'capped at ' + str(cap) if capped else 'exact'})")

        converged = set(new_medoids) == set(medoids)
        medoids = new_medoids
        clusters = new_clusters
        if converged:
            print(f"  Clustering converged after {iteration + 1} iteration(s).")
            break
    else:
        print(f"  Clustering reached max_iter={max_iter} without convergence.")

    return [clusters[m] for m in medoids if clusters.get(m)]


# ═══════════════════════════════════════════════════════════════
# Binned reverse translation
# ═══════════════════════════════════════════════════════════════

def run_binned(seqs, k, sample_size, seed, pick_fn, verbose=False,
               overall_pb=None, overall_done=None, overall_total=None,
               verbose_fh=None, medoid_cap=2000):
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
    bins = cluster_sequences(working, k=k, seed=seed, medoid_cap=medoid_cap)
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
# Stop-risk sweep  (w_stop tradeoff table)
# ═══════════════════════════════════════════════════════════════

def run_stop_sweep(seqs, args, verbose_fh=None):
    """
    Sweep w_stop over a range at fixed w_loss/w_offtarget, clustering
    once (clustering doesn't depend on w_stop) and re-running codon
    selection per bin at each w_stop value via best_codon_stop_aware
    over the unrestricted CODON_AA_MAP_ALL.

    Prints one table per bin: avg loss, avg off-target, and the true
    sequence-level stop-free risk (from stop_risk_estimate) at each
    w_stop -- a Pareto tradeoff curve to read a value off of, rather
    than a single point estimate you'd otherwise have to guess at.

    verbose_fh: if provided and args.verbose is set, per-position detail
                for every (bin, w_stop) combination is written there
                with a section header, same pattern as run_binned.

    Returns csv_rows: a list of per-(bin, w_stop) dicts, for --csv output.
    """
    k = args.bins
    n_orig = len(seqs)

    if args.sample_size > 0 and args.sample_size < n_orig:
        rng = random.Random(args.seed)
        working = rng.sample(seqs, args.sample_size)
        print(f"  Sampled {args.sample_size} / {n_orig} sequences for clustering.")
    else:
        working = seqs
        print(f"  Using all {n_orig} sequences (no sampling).")

    if k == 1:
        bins = [working]
    else:
        print(f"  Clustering into k={k} bins...")
        bins = cluster_sequences(working, k=k, seed=args.seed, medoid_cap=args.medoid_cap)

    wl  = args.w_loss if args.w_loss is not None else 1.0
    wot = args.w_offtarget if args.w_offtarget is not None else 1.0
    sweep_ws = _make_weight_range(args.stop_sweep_max, args.stop_sweep_step)

    print(SEP)
    print(f"STOP-RISK SWEEP  (w_loss={wl}  w_offtarget={wot}  "
          f"w_stop = 0..{args.stop_sweep_max} step {args.stop_sweep_step})")
    print("  Searching CODON_AA_MAP_ALL (unrestricted) at every w_stop value.")
    print(SEP + "\n")

    csv_rows = []
    actual_k = len(bins)
    for i, bin_seqs in enumerate(bins, 1):
        cols = columns(bin_seqs)
        print(f"  Bin {i}/{actual_k}  ({len(bin_seqs)} sequences)")
        print(f"  {'w_stop':>8}  {'avg loss':>10}  {'avg off-tgt':>12}  {'seq stop risk':>14}")
        print(f"  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*14}")
        for w_stop in sweep_ws:
            pick = lambda s, wl=wl, wot=wot, ws=w_stop: best_codon_stop_aware(
                s, wl, wot, ws, CODON_AA_MAP_ALL)

            vfile_to_use = None
            if args.verbose and verbose_fh is not None:
                vfile_to_use = verbose_fh
                verbose_fh.write(f"\n{'=' * 60}\n")
                verbose_fh.write(f"BIN {i}/{actual_k}  w_stop={w_stop:.2f}"
                                 f"  ({len(bin_seqs)} sequences)\n")
                verbose_fh.write(f"{'=' * 60}\n\n")

            results, sequence = run_pass(
                cols, pick, verbose=args.verbose,
                label=f"bin{i} w_stop={w_stop:.2f}",
                verbose_file=vfile_to_use,
                show_progress=False,
            )
            avg_loss = sum(r["loss"] for r in results) / len(results) if results else 0.0
            avg_ot   = sum(r["off_target"] for r in results) / len(results) if results else 0.0
            p_stop, _ = stop_risk_estimate(results)
            print(f"  {w_stop:>8.2f}  {avg_loss:>10.4f}  {avg_ot:>12.4f}  {p_stop:>14.4f}")

            csv_rows.append({
                "w_loss":              wl,
                "w_offtarget":         wot,
                "w_stop":              w_stop,
                "bin":                 i,
                "bin_size":            len(bin_seqs),
                "avg_loss":            round(avg_loss, 6),
                "avg_offtarget":       round(avg_ot, 6),
                "seq_stop_risk":       round(p_stop, 6),
                "degenerate_sequence": sequence,
            })
        print()

    print("  Interpretation:")
    print("  - w_stop=0 is the plain w_loss/w_offtarget optimum over the unrestricted map.")
    print("  - 'seq stop risk' is the true P(>=1 spontaneous stop) for that bin's oligo,")
    print("    i.e. the same number stop_risk_estimate() reports elsewhere.")
    print("  - Pick the smallest w_stop that gets seq stop risk where you need it;")
    print("    pushing w_stop further than that only costs recall/precision for no benefit.")
    print()

    return csv_rows


# ═══════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════

SEP = "=" * 60

def parse_args():
    p = argparse.ArgumentParser(
        description="Binned degenerate IUPAC codon design, stop-codon-safe by default (V5).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("fasta",
        help="Aligned FASTA file of protein sequences")

    # ── Stop-codon handling (V5) ────────────────────────────────
    p.add_argument("--allow-stop-codons", action="store_true",
        help="Search the full unrestricted IUPAC codon space (V4 "
             "behaviour), which may spontaneously encode a stop codon "
             "for some fraction of synthesized oligos. Default: OFF — "
             "only codons whose expansion is 100%% stop-free are used.")

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
    p.add_argument("--medoid-cap", type=int, default=2000, metavar="N",
        help="CLARA-style cap on cluster size used when recomputing a "
             "medoid each clustering iteration (default: 2000). Bounds "
             "the O(cluster_size^2) update step so --sample-size 0 stays "
             "tractable on large inputs. Set 0 to disable (exact, but "
             "slow on very large clusters).")
    p.add_argument("--n-jobs", type=int, default=1, metavar="N",
        help="Worker processes for parallel grid search in --mode grid "
             "/ --mode surface (default: 1 = sequential). Each weight "
             "combination is independent, so this parallelizes cleanly; "
             "has no effect in --mode direct or --bins-sweep.")

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
    p.add_argument("--w-stop", type=float, default=None, metavar="F",
        help="Manual stop-codon-risk weight (supply together with "
             "--w-loss and --w-offtarget). Searches the unrestricted "
             "codon space (CODON_AA_MAP_ALL) and penalizes each "
             "candidate by -log(1 - stop_frac), so the position-by-"
             "position sum exactly equals -log(true sequence-level "
             "stop-free probability) rather than a linear approximation "
             "of it. w_stop=0 behaves like plain --w-loss/--w-offtarget "
             "over the unrestricted map.")
    p.add_argument("--stop-sweep", action="store_true",
        help="Instead of a single run, sweep w_stop over a range (see "
             "--stop-sweep-max/--stop-sweep-step) at fixed "
             "w_loss/w_offtarget (from --w-loss/--w-offtarget, default "
             "1.0/1.0 if not given) and print a table per bin of avg "
             "loss / avg off-target / true sequence-level stop risk at "
             "each w_stop -- a tradeoff curve rather than one answer. "
             "Overrides --mode and --w-stop.")
    p.add_argument("--stop-sweep-max", type=float, default=2.0, metavar="F",
        help="Upper bound of the w_stop sweep range for --stop-sweep (default: 2.0)")
    p.add_argument("--stop-sweep-step", type=float, default=0.2, metavar="F",
        help="Step size of the w_stop sweep range for --stop-sweep (default: 0.2)")

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
    if args.w_stop is not None and (args.w_loss is None or args.w_offtarget is None):
        p.error("--w-stop requires --w-loss and --w-offtarget to also be supplied.")
    if args.w_stop is not None and args.w_stop < 0:
        p.error("--w-stop must be >= 0.")
    if args.stop_sweep and args.w_stop is not None:
        p.error("--stop-sweep sweeps w_stop itself; do not also supply --w-stop.")
    if args.stop_sweep and args.bins_sweep:
        p.error("--stop-sweep and --bins-sweep cannot both be set.")
    if args.stop_sweep and args.sweep_max != 5:
        p.error("--sweep-max only applies to --bins-sweep and has no effect on "
                 "--stop-sweep (which sweeps w_stop, not k). Did you mean "
                 "--stop-sweep-max?")
    if args.stop_sweep and args.stop_sweep_step <= 0:
        p.error("--stop-sweep-step must be > 0.")
    if args.stop_sweep and args.stop_sweep_max < 0:
        p.error("--stop-sweep-max must be >= 0.")

    return args


def main():
    args = parse_args()

    # ── Start log tee before any output ────────────────────────
    tee = None
    if args.log:
        tee = Tee(args.log)
        # Write a header so the log is self-contained
        import datetime
        print(f"# reverse_translation_V5_StopSafe.py  log")
        print(f"# Date : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"# Args : {' '.join(sys.argv[1:])}")
        print(f"# {'=' * 56}")

    seqs = read_fasta(args.fasta)
    print(f"\nLoaded {len(seqs)} sequences, length {len(seqs[0])} AA each.")

    # ── Resolve active codon search space (V5) ──────────────────
    if args.allow_stop_codons:
        active_codon_map = CODON_AA_MAP_ALL
        codon_mode_str    = ("unsafe — all IUPAC codons (V4 behaviour); "
                              "some may spontaneously form a stop codon")
    else:
        active_codon_map = CODON_AA_MAP_SAFE
        codon_mode_str    = "stop-safe (default) — only codons with zero stop-codon risk"
    print(f"Codon search space : {codon_mode_str}")
    print(f"  ({len(active_codon_map)} degenerate codons available"
          f"{'' if args.allow_stop_codons else f', vs {len(CODON_AA_MAP_ALL)} if --allow-stop-codons were set'})")

    # ── Resolve pick function ───────────────────────────────────
    # Manual weights override --mode entirely (same as V3).
    # In sweep mode, always use direct regardless of --mode.
    if args.w_stop is not None:
        wl, wot, ws = args.w_loss, args.w_offtarget, args.w_stop
        # --w-stop always searches the unrestricted map: every codon in
        # CODON_AA_MAP_SAFE already has zero stop risk, so w_stop would
        # have nothing to trade against if we stayed on active_codon_map
        # when --allow-stop-codons wasn't also set.
        pick_fn  = lambda s: best_codon_stop_aware(s, wl, wot, ws, CODON_AA_MAP_ALL)
        mode_str = f"manual weights (stop-aware)  w_loss={wl}  w_offtarget={wot}  w_stop={ws}"
        if not args.allow_stop_codons:
            print("  Note: --w-stop searches the unrestricted codon space "
                  "(CODON_AA_MAP_ALL) regardless of --allow-stop-codons, "
                  "since w_stop needs stop-risky codons available to "
                  "actually trade against.")
    elif args.w_loss is not None:
        wl, wot = args.w_loss, args.w_offtarget
        pick_fn  = lambda s: best_codon_weighted(s, wl, wot, active_codon_map)
        mode_str = f"manual weights  w_loss={wl}  w_offtarget={wot}"
    else:
        pick_fn  = lambda s: best_codon_direct(s, active_codon_map)
        mode_str = "direct"   # overridden below for grid/surface in single-run path

    # ── Open single verbose output file if requested ──────────
    verbose_fh = None
    if args.verbose:
        import os
        fasta_base = os.path.splitext(os.path.basename(args.fasta))[0]
        verbose_path = f"verbose_{fasta_base}.txt"
        verbose_fh = open(verbose_path, "w", encoding="utf-8")
        verbose_fh.write(f"Verbose per-position detail\n")
        verbose_fh.write(f"Input : {args.fasta}\n")
        verbose_fh.write(f"Mode  : {'stop-sweep' if args.stop_sweep else ('bins-sweep' if args.bins_sweep else args.mode)}\n")
        verbose_fh.write(f"Bins  : {'sweep k=1..5' if args.bins_sweep else args.bins}\n")
        verbose_fh.write("=" * 60 + "\n")
        print(f"  Verbose detail will be written to: {verbose_path}")

    if args.stop_sweep:
        print(SEP)
        print(f"MODE: stop-risk sweep  —  k = {args.bins}")
        print(f"Sample size : {args.sample_size if args.sample_size > 0 else 'all'}")
        print(SEP + "\n")
        csv_rows = run_stop_sweep(seqs, args, verbose_fh=verbose_fh)

    elif args.bins_sweep:
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
                pick_fn=lambda s: best_codon_direct(s, active_codon_map),
                verbose=args.verbose,
                overall_pb=overall_pb,
                overall_done=overall_done,
                overall_total=total_passes,
                verbose_fh=verbose_fh if args.verbose and k > 1 else None,
                medoid_cap=args.medoid_cap,
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
                stop_stat_line = stop_codon_summary_line(
                    columns(bin_seqs), results, args.allow_stop_codons
                )
                print(f"  Bin {i} ({bin_size} sequences):")
                print(f"    Degenerate sequence: {sequence}")
                print(f"    Encoded space      : {space_str}")
                print(f"    Seqs covered       : {n_cov}/{n_tot}  ({frac*100:.2f}%)")
                print(f"    Stop-codon stat    : {stop_stat_line}")
                csv_rows.append({
                    "k":                    k,
                    "bin":                  i,
                    "bin_size":             bin_size,
                    "encoded_space":        bin_space,
                    "pct_of_baseline":      round(pct, 4),
                    "reduction_pct":        round(reduction_bin, 4),
                    "seqs_covered_pct":     round(frac * 100, 4),
                    "degenerate_sequence":  sequence,
                    "stop_codon_stat":      stop_stat_line,
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
            pick_fn=lambda s: best_codon_direct(s, active_codon_map),
            verbose=False,
            medoid_cap=args.medoid_cap,
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
                bins = cluster_sequences(working, k=k, seed=args.seed, medoid_cap=args.medoid_cap)

            bin_results = []
            total_space_k = 0

            for i, bin_seqs in enumerate(bins, 1):
                actual_k = len(bins)
                print(f"\n{SEP}")
                print(f"  Bin {i}/{actual_k}  ({len(bin_seqs)} sequences) — {args.mode} mode")
                print(SEP)
                cols = columns(bin_seqs)

                # Direct baseline for this bin
                res_direct, seq_direct = run_pass(
                    cols, lambda s: best_codon_direct(s, active_codon_map), verbose=False,
                    label=f"bin{i} direct  ")
                direct_obj = sum(r["direct_obj"] for r in res_direct) / len(res_direct)
                print(f"  Direct optimum  avg(2−R−P) = {direct_obj:.6f}\n")

                n_combos = round((args.grid_max / args.grid_step + 1) ** 2)
                print(f"  Running {n_combos} weight combinations...")
                grid_res, best_w, best_grid_obj = grid_search(
                    cols, args.grid_max, args.grid_step, active_codon_map, n_jobs=args.n_jobs)
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
                    pick = lambda s, wl=best_w[0], wot=best_w[1]: best_codon_weighted(
                        s, wl, wot, active_codon_map)
                    results, sequence = run_pass(cols, pick, verbose=args.verbose,
                                                  label=f"bin{i} grid    ",
                                                  verbose_file=vfile_to_use)
                    summarise(results, label=f"bin {i} (grid best weights)", sequence=sequence,
                              baseline_space=total_space_baseline, coverage_seqs=working)
                    print_stop_codon_stats(cols, results, args.allow_stop_codons, label=f"bin {i}")

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
                    pick  = lambda s, wl=x_use, wot=y_use: best_codon_weighted(
                        s, wl, wot, active_codon_map)
                    results, sequence = run_pass(cols, pick, verbose=args.verbose,
                                                  label=f"bin{i} surface ",
                                                  verbose_file=vfile_to_use)
                    summarise(results, label=f"bin {i} (surface weights)", sequence=sequence,
                              baseline_space=total_space_baseline, coverage_seqs=working)
                    print_stop_codon_stats(cols, results, args.allow_stop_codons, label=f"bin {i}")

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
                medoid_cap=args.medoid_cap,
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
                print_stop_codon_stats(columns(bin_seqs), results, args.allow_stop_codons,
                                        label=f"bin {i}")

            # Also show direct comparison if manual weights were used
            if args.w_loss is not None:
                res_direct, _ = run_pass(
                    columns(seqs[:min(len(seqs), args.sample_size or len(seqs))]),
                    lambda s: best_codon_direct(s, active_codon_map),
                    verbose=False, label="direct check  ")
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
        print(f"  Codon search space                   : {codon_mode_str}")
        print(f"  k=1 encoded space (direct baseline) : {total_space_baseline}")
        print(f"  k={k} total encoded space            : {total_space_k}")
        print(f"  Space reduction                      : {reduction:.2f}%")
        print(f"  Oligos to order                      : {len(bin_results)}  (vs 1 for k=1)")
        print()


    # ── Write CSV if requested ─────────────────────────────────
    if args.csv and args.stop_sweep:
        import csv
        fieldnames = ["w_loss", "w_offtarget", "w_stop", "bin", "bin_size",
                      "avg_loss", "avg_offtarget", "seq_stop_risk",
                      "degenerate_sequence"]
        with open(args.csv, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"  CSV written to: {args.csv}")
    elif args.csv and args.bins_sweep:
        import csv
        fieldnames = ["k", "bin", "bin_size", "encoded_space",
                      "pct_of_baseline", "reduction_pct", "seqs_covered_pct",
                      "degenerate_sequence", "stop_codon_stat"]
        with open(args.csv, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"  CSV written to: {args.csv}")
    elif args.csv:
        print("  Note: --csv is only populated during --bins-sweep or --stop-sweep runs.")

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