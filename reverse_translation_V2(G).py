"""
bio_model_v2.py  –  degenerate codon design with optimised loss & off-target metrics

Usage
-----
    python bio_model_v2.py <alignment.fasta> [--w-loss W] [--w-offtarget W]

    --w-loss       Weight for recall loss   (default 1.0)
    --w-offtarget  Weight for off-target rate (default 0.5)

    Increasing --w-offtarget pushes the optimiser toward cleaner codons (fewer
    junk AAs) at the cost of potentially missing more of the observed diversity.
    Increasing --w-loss does the opposite.
"""

from collections import Counter
from itertools import product
import sys
import math
import argparse

# =========================================================
# 3-group biochemical model
# =========================================================

GROUPS = {
    "hydrophobic": {"A","V","L","I","M","F","W","Y","P"},
    "polar":       {"S","T","N","Q","C","G"},
    "charged":     {"D","E","K","R","H"},
}

AA_TO_GROUP = {
    aa: g
    for g, aas in GROUPS.items()
    for aa in aas
}

# =========================================================
# Standard genetic code (full 64-codon table)
# =========================================================

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

# =========================================================
# IUPAC expansion
# =========================================================

IUPAC = {
    "A":{"A"},"C":{"C"},"G":{"G"},"T":{"T"},
    "R":{"A","G"},"Y":{"C","T"},"S":{"G","C"},
    "W":{"A","T"},"K":{"G","T"},"M":{"A","C"},
    "B":{"C","G","T"},"D":{"A","G","T"},
    "H":{"A","C","T"},"V":{"A","C","G"},"N":{"A","C","G","T"},
}

def expand(codon):
    """Expand an IUPAC codon to all concrete DNA triplets it represents."""
    return map("".join, product(*[IUPAC[b] for b in codon]))

# =========================================================
# Pre-compute codon → frozenset of amino acids (no stops)
# Doing this once at import is ~100× faster than re-expanding every call
# =========================================================

def _build_codon_aa_map():
    mapping = {}
    for c in ("".join(p) for p in product(IUPAC.keys(), repeat=3)):
        aas = frozenset(
            CODON_TABLE[x]
            for x in expand(c)
            if x in CODON_TABLE and CODON_TABLE[x] != "*"
        )
        if aas:                      # skip all-stop codons
            mapping[c] = aas
    return mapping

CODON_AA_MAP = _build_codon_aa_map()   # {iupac_codon: frozenset of AAs}

# =========================================================
# FASTA parser
# =========================================================

def read_fasta(file):
    seqs, cur = [], []

    with open(file) as f:
        for line in f:
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

# =========================================================
# Column extraction
# =========================================================

def columns(seqs):
    """Per-column amino acid lists, ignoring gap characters."""
    L = len(seqs[0])
    return [
        [s[i] for s in seqs if s[i] != "-"]
        for i in range(L)
    ]

# =========================================================
# Shannon entropy (bits)
# =========================================================

def entropy(col):
    c = Counter(col)
    total = len(col)
    return -sum(
        (v / total) * math.log2(v / total)
        for v in c.values() if v > 0
    )

# =========================================================
# Per-position metrics
# =========================================================

def position_metrics(orig_set, enc_aas):
    """
    orig_set  : set of AAs observed in the alignment column
    enc_aas   : frozenset of AAs encoded by the chosen IUPAC codon

    Returns
    -------
    recall        : fraction of observed AAs that are captured  (= 1 - loss)
    precision     : fraction of encoded AAs that were actually observed
    loss          : 1 - recall   (missed targets)
    off_target    : 1 - precision  (junk / unwanted AAs encoded)
    f1            : harmonic mean of precision and recall
    """
    overlap   = orig_set & enc_aas
    recall    = len(overlap) / len(orig_set)           if orig_set  else 1.0
    precision = len(overlap) / len(enc_aas)            if enc_aas   else 0.0
    loss      = 1.0 - recall
    off_target = 1.0 - precision
    denom     = precision + recall
    f1        = 2 * precision * recall / denom         if denom > 0 else 0.0
    return recall, precision, loss, off_target, f1

# =========================================================
# Optimised codon search
#
# Objective (minimise):
#   w_loss * loss  +  w_offtarget * off_target_rate
#
# This is exhaustive over all 3,375 IUPAC codons so the result is
# globally optimal for the chosen weights.
# =========================================================

def best_codon(orig_set, w_loss=1.0, w_offtarget=0.5):
    """
    Find the IUPAC codon that minimises the weighted objective.

    Parameters
    ----------
    orig_set     : set of AAs observed at this alignment column
    w_loss       : penalty weight for missing observed AAs (recall loss)
    w_offtarget  : penalty weight for encoding unwanted AAs

    Returns
    -------
    best_codon   : str
    best_enc_aas : frozenset
    best_obj     : float  (objective value)
    """
    best_codon_str = None
    best_enc_aas   = frozenset()
    best_obj       = float("inf")

    for codon, enc_aas in CODON_AA_MAP.items():
        _, _, loss, off_target, _ = position_metrics(orig_set, enc_aas)
        obj = w_loss * loss + w_offtarget * off_target

        if obj < best_obj:
            best_obj       = obj
            best_codon_str = codon
            best_enc_aas   = enc_aas

    return best_codon_str, best_enc_aas, best_obj

# =========================================================
# MAIN
# =========================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Design a degenerate IUPAC codon sequence from a protein MSA."
    )
    p.add_argument("fasta", help="Aligned FASTA file (protein sequences)")
    p.add_argument(
        "--w-loss", type=float, default=1.0, metavar="W",
        help="Penalty weight for recall loss — missing observed AAs (default: 1.0)"
    )
    p.add_argument(
        "--w-offtarget", type=float, default=0.5, metavar="W",
        help="Penalty weight for off-target AAs encoded (default: 0.5)"
    )
    return p.parse_args()


def main():
    args = parse_args()

    print(f"Weights  →  loss: {args.w_loss}  |  off-target: {args.w_offtarget}\n")

    seqs = read_fasta(args.fasta)
    cols = columns(seqs)

    # Accumulators
    original_space   = 1
    encoded_space    = 1
    off_target_aas   = 0    # total unwanted AAs summed across positions
    total_encoded    = 0    # total AAs encoded across positions (denominator)
    metrics_per_pos  = []
    final            = []

    for i, col in enumerate(cols):

        orig_set = set(col) - {"-"}
        if not orig_set:
            print(f"Pos {i+1}: all gaps — skipping\n")
            final.append("NNN")
            continue

        H = entropy(col)

        codon, enc_aas, obj = best_codon(
            orig_set,
            w_loss=args.w_loss,
            w_offtarget=args.w_offtarget,
        )

        recall, precision, loss, off_target, f1 = position_metrics(orig_set, enc_aas)

        # Space metrics
        original_space *= max(len(orig_set), 1)
        encoded_space  *= max(len(enc_aas),  1)

        # Off-target AA counts (for global off-target efficiency below)
        junk_aas      = enc_aas - orig_set
        off_target_aas += len(junk_aas)
        total_encoded  += len(enc_aas)

        metrics_per_pos.append((recall, precision, loss, off_target, f1))
        final.append(codon)

        print(f"Pos {i+1}")
        print(f"  Entropy (bits)   : {H:.3f}")
        print(f"  Observed AAs     : {sorted(orig_set)}")
        print(f"  Encoded AAs      : {sorted(enc_aas)}")
        print(f"  Off-target AAs   : {sorted(junk_aas)}")
        print(f"  Chosen codon     : {codon}")
        print(f"  Recall (1-loss)  : {recall:.3f}")
        print(f"  Precision        : {precision:.3f}")
        print(f"  Off-target rate  : {off_target:.3f}")
        print(f"  F1 score         : {f1:.3f}")
        print(f"  Objective value  : {obj:.4f}\n")

    if not metrics_per_pos:
        print("No valid positions found.")
        return

    # ── Global metrics ─────────────────────────────────────────────────────

    # Compression: how much smaller is the encoded space vs. the original?
    compression = (
        original_space / encoded_space if encoded_space > 0 else float("inf")
    )

    # Off-target efficiency: what fraction of all encoded AAs across the
    # full library are unwanted?  0 = perfect (no junk), 1 = everything junk
    global_off_target = (
        off_target_aas / total_encoded if total_encoded > 0 else 0.0
    )

    avg = lambda idx: sum(m[idx] for m in metrics_per_pos) / len(metrics_per_pos)
    avg_recall     = avg(0)
    avg_precision  = avg(1)
    avg_loss       = avg(2)
    avg_off_target = avg(3)
    avg_f1         = avg(4)

    print("=" * 50)
    print("FINAL RESULTS")
    print("=" * 50)

    print("\nDegenerate sequence:")
    print("".join(final))

    print("\n── Per-position averages ──────────────────────")
    print(f"  Avg recall (coverage)    : {avg_recall:.4f}")
    print(f"  Avg loss (missed targets): {avg_loss:.4f}")
    print(f"  Avg precision            : {avg_precision:.4f}")
    print(f"  Avg off-target rate      : {avg_off_target:.4f}")
    print(f"  Avg F1                   : {avg_f1:.4f}")

    print("\n── Library-level metrics ──────────────────────")
    print(f"  Compression ratio        : {compression:.4f}")
    print(f"    (original AA space / encoded AA space)")
    print(f"  Global off-target rate   : {global_off_target:.4f}")
    print(f"    (unwanted AAs / total encoded AAs across all positions)")

    print("\n── Interpretation ─────────────────────────────")
    print("  Recall → 1.0   : all observed AAs are covered")
    print("  Loss   → 0.0   : no observed AAs are missed")
    print("  Precision → 1.0: no junk AAs encoded")
    print("  Off-target → 0.0: library is clean")
    print("  Compression > 1: library is smaller than naive enumeration")
    print("  Adjust --w-loss / --w-offtarget to shift the tradeoff.\n")


if __name__ == "__main__":
    main()