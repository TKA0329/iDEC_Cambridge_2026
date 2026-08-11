#!/usr/bin/env python3
"""
IDP Solubility Tag Sequence Generator
======================================
In silico de novo generation of IDP solubility tag candidate sequences.

Pipeline:
  Stage 1 — Motif Generation  (KPEV+, Motif2, 5AA, 10AA, Hydrophobic Patches)
  Stage 2 — Repeat Algorithm
  Stage 3 — Recombination Algorithm

Output: CSV with sequences and full provenance metadata.
"""

import csv
import random
import itertools
import os
import sys
import time
import hashlib
from typing import List, Dict, Tuple, Generator, Optional

# ================================================================
#  CONFIGURATION — All tunable parameters live here
# ================================================================
CONFIG = {

    # ── Samples per motif pattern ──────────────────────────────
    'sample_counts': {
        'motif2':        10,   # Motif 2  (4 AA)
        '5aa':           10,   # 5 AA motifs
        '10aa':           3,   # 10 AA motifs (reduced to manage scale)
        'recombination': 10,   # combinations sampled per slot-config × length
    },

    # ── Amino acid ratios  {AA: weight} ───────────────────────
    #    Change weights to shift sampling probability.
    #    Add/remove AAs to change the alphabet.
    'ratios': {
        'positive': {'K': 1, 'R': 1},
        'negative': {'E': 1, 'D': 1},
        'np1':      {'G': 1, 'P': 3},          # non-polar set 1 (motifs 3 & 4)
        'np2':      {'G': 2, 'A': 1, 'V': 1,   # non-polar set 2 (motif 2)
                     'I': 1, 'L': 1},
        'polar':    {'S': 2, 'T': 2, 'Q': 2,   # includes H (conditionally charged)
                     'N': 1, 'H': 2},
    },

    # ── Target lengths and tolerance ──────────────────────────
    'target_lengths':   [15, 50, 100, 150, 200, 250],
    'length_tolerance': 3,

    # ── Stochastic set parameters ─────────────────────────────
    'set_attempts':      3,     # repeats for Set 2 and Set 3
    'set3_remove_prob':  0.75,  # probability of removing a +  (Set 3)
    # replace probability is implicitly 1 - set3_remove_prob

    # ── I/O ───────────────────────────────────────────────────
    'output_file':      'generated_sequences.csv',
    'write_batch_size': 10_000,   # rows buffered before flushing to disk

    # ── Reproducibility ───────────────────────────────────────
    'random_seed': 42,   # set to None for a different run each time
}

# ================================================================
#  WEIGHTED SAMPLING
# ================================================================

_POOLS: Dict[str, List[str]] = {}


def _build_pool(ratio_dict: Dict[str, int]) -> List[str]:
    pool = []
    for aa, weight in ratio_dict.items():
        pool.extend([aa] * weight)
    return pool


def _init_pools() -> None:
    global _POOLS
    _POOLS = {cat: _build_pool(r) for cat, r in CONFIG['ratios'].items()}


def _pick(category: str) -> str:
    """Sample one amino acid from a category using configured weights."""
    return random.choice(_POOLS[category])


# ================================================================
#  STAGE 1 — MOTIF GENERATION
# ================================================================

def gen_kpev() -> List[str]:
    """Motif 1 — Fixed biological anchor (Smt3-based). No sampling."""
    return ['KPEV+']


def gen_motif2(n: int) -> List[str]:
    """
    Motif 2 — (pos)–P–(neg)–np2–+
    Pattern:  (K|R) – P – (E|D) – (G|A|V|I|L) – +
    Samples n unique sequences using configured ratios.
    """
    seen: set = set()
    max_tries = n * 500
    tries = 0
    while len(seen) < n and tries < max_tries:
        m = _pick('positive') + 'P' + _pick('negative') + _pick('np2') + '+'
        seen.add(m)
        tries += 1
    if len(seen) < n:
        print(f"  ⚠  Motif 2: requested {n} unique sequences, "
              f"achieved {len(seen)} (alphabet too small).")
    return list(seen)


def gen_5aa_motifs(n: int) -> List[str]:
    """
    5 AA motifs — ends with +.
    Rules:
      • Exactly 1 np1 slot  (position varies across all 5 positions)
      • Remaining 4 slots each independently from {neg, pos, polar}
      • n samples taken per (np_position × category_combo) pattern
    Total patterns: 5 × 3^4 = 405
    """
    cats = ['negative', 'positive', 'polar']
    seen: set = set()

    for np_pos in range(5):
        for combo in itertools.product(cats, repeat=4):
            for _ in range(n):
                seq = []
                ci = 0
                for pos in range(5):
                    if pos == np_pos:
                        seq.append(_pick('np1'))
                    else:
                        seq.append(_pick(combo[ci]))
                        ci += 1
                seen.add(''.join(seq) + '+')
    return list(seen)


def gen_10aa_motifs(n: int) -> List[str]:
    """
    10 AA motifs — ends with +.
    Rules:
      • Exactly 1 np1 slot  (position varies across all 10 positions)
      • Remaining 9 slots each independently from {neg, pos, polar}
      • n samples taken per (np_position × category_combo) pattern
    Total patterns: 10 × 3^9 = 196,830
    """
    cats = ['negative', 'positive', 'polar']
    seen: set = set()
    total_patterns = 10 * (3 ** 9)  # 196,830
    done = 0
    t0 = time.time()

    for np_pos in range(10):
        for combo in itertools.product(cats, repeat=9):
            for _ in range(n):
                seq = []
                ci = 0
                for pos in range(10):
                    if pos == np_pos:
                        seq.append(_pick('np1'))
                    else:
                        seq.append(_pick(combo[ci]))
                        ci += 1
                seen.add(''.join(seq) + '+')
            done += 1
            if done % 20_000 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (total_patterns - done) / rate if rate > 0 else 0
                print(
                    f"\r    {done:>7,}/{total_patterns:,} patterns "
                    f"({100*done/total_patterns:.1f}%)  "
                    f"ETA: {eta:.0f}s   ",
                    end='', flush=True,
                )
    print()
    return list(seen)


def gen_hydrophobic_patches() -> List[str]:
    """
    All 25 ordered 2 AA combinations of {G, A, V, I, L}.
    No '+' — used as patch replacements at + positions.
    """
    aa = ['G', 'A', 'V', 'I', 'L']
    return [a + b for a in aa for b in aa]


# ================================================================
#  SET PROCESSING  (Set 1 / Set 2 / Set 3)
# ================================================================

def apply_set1(raw: str) -> str:
    """Set 1 — Remove all '+' symbols. Deterministic."""
    return raw.replace('+', '')


def apply_set2(raw: str, patches: List[str]) -> str:
    """Set 2 — Replace every '+' with a randomly chosen hydrophobic patch."""
    return ''.join(random.choice(patches) if c == '+' else c for c in raw)


def apply_set3(raw: str, patches: List[str], remove_prob: float) -> str:
    """
    Set 3 — Each '+' is independently:
      • Removed   with probability remove_prob
      • Replaced  with a random patch with probability (1 - remove_prob)
    """
    out = []
    for c in raw:
        if c == '+':
            if random.random() < remove_prob:
                pass                            # remove
            else:
                out.append(random.choice(patches))  # replace
        else:
            out.append(c)
    return ''.join(out)


def triplication(
    raw: str,
    patches: List[str],
    attempts: int,
    remove_prob: float,
) -> List[Tuple[str, int, int]]:
    """
    Apply the three-set triplication rule.
    Returns list of (resolved_sequence, set_number, attempt_number).

    Set 1 — 1 result  (deterministic, not triplicated)
    Set 2 — `attempts` results
    Set 3 — `attempts` results
    """
    results = [(apply_set1(raw), 1, 1)]
    for i in range(1, attempts + 1):
        results.append((apply_set2(raw, patches), 2, i))
    for i in range(1, attempts + 1):
        results.append((apply_set3(raw, patches, remove_prob), 3, i))
    return results


# ================================================================
#  LENGTH HELPERS
# ================================================================

def valid_repeat_counts(aa_len: int, target: int, tol: int) -> List[int]:
    """
    Return all n ≥ 1 such that  (target - tol) ≤ n × aa_len ≤ (target + tol).
    """
    lo, hi = target - tol, target + tol
    if aa_len == 0:
        return []
    return [n for n in range(1, hi // aa_len + 2) if lo <= n * aa_len <= hi]


def valid_mixed_slot_configs(
    targets: List[int], tol: int
) -> Dict[int, List[Tuple[int, int]]]:
    """
    For each target length, find all (a, b) with a > 0 and b > 0 such that
    a × 5 + b × 10  is within [target - tol, target + tol].
    Pure-5AA (b=0) and pure-10AA (a=0) are excluded (handled separately).
    """
    result: Dict[int, List[Tuple[int, int]]] = {}
    for t in targets:
        lo, hi = t - tol, t + tol
        combos = []
        for b in range(1, hi // 10 + 1):
            for a in range(1, (hi - b * 10) // 5 + 1):
                if lo <= a * 5 + b * 10 <= hi:
                    combos.append((a, b))
        result[t] = combos
    return result


# ================================================================
#  STAGE 2 — REPEAT ALGORITHM
# ================================================================

def repeat_algorithm(
    motifs:      List[str],
    motif_type:  str,
    targets:     List[int],
    tol:         int,
    patches:     List[str],
    attempts:    int,
    remove_prob: float,
) -> Generator[Dict, None, None]:
    """
    For every motif × valid target length:
      1. Repeat the motif n times to reach the target length.
      2. Apply Set 1 / 2 / 3 triplication.
    Yields one provenance dict per output sequence.
    """
    for motif in motifs:
        aa_len = len(motif) - 1          # exclude trailing '+'
        for target in targets:
            for n in valid_repeat_counts(aa_len, target, tol):
                raw = motif * n           # e.g. "KPEV+KPEV+KPEV+"
                for seq, set_num, attempt in triplication(
                    raw, patches, attempts, remove_prob
                ):
                    yield {
                        'source':        'repeat',
                        'motif_type':    motif_type,
                        'motif_pattern': motif,
                        'n_repeats':     n,
                        'target_length': target,
                        'set':           set_num,
                        'attempt':       attempt,
                        'sequence':      seq,
                        'actual_length': len(seq),
                    }


# ================================================================
#  STAGE 3 — RECOMBINATION ALGORITHM
# ================================================================

def _sample_combination(pool: List[str], n_slots: int) -> str:
    """Pick n_slots motifs uniformly from pool and concatenate (with embedded '+')."""
    return ''.join(random.choice(pool) for _ in range(n_slots))


def recombination_algorithm(
    pool_4aa:    List[str],
    pool_5aa:    List[str],
    pool_10aa:   List[str],
    targets:     List[int],
    tol:         int,
    patches:     List[str],
    attempts:    int,
    remove_prob: float,
    n_samples:   int,
) -> Generator[Dict, None, None]:
    """
    Four recombination categories per target length:
      4AA-only  — Motif 2 pool recombined with itself
      5AA-only  — 5AA pool recombined with itself
      10AA-only — 10AA pool recombined with itself
      mixed     — 5AA + 10AA interleaved (all valid (a,b) slot configs)

    KPEV+ (Motif 1) is excluded from recombination.

    For each category × target length × slot config:
      → sample n_samples combinations
      → apply Set 1 / 2 / 3 triplication to each
    """
    mixed_configs = valid_mixed_slot_configs(targets, tol)

    for target in targets:

        # ── 4AA recombination ──────────────────────────────────
        if pool_4aa:
            for n in valid_repeat_counts(4, target, tol):
                for _ in range(n_samples):
                    raw = _sample_combination(pool_4aa, n)
                    for seq, set_num, attempt in triplication(
                        raw, patches, attempts, remove_prob
                    ):
                        yield {
                            'source':        'recombination',
                            'motif_type':    '4aa',
                            'motif_pattern': raw,
                            'n_repeats':     n,
                            'target_length': target,
                            'set':           set_num,
                            'attempt':       attempt,
                            'sequence':      seq,
                            'actual_length': len(seq),
                        }

        # ── 5AA recombination ──────────────────────────────────
        if pool_5aa:
            for n in valid_repeat_counts(5, target, tol):
                for _ in range(n_samples):
                    raw = _sample_combination(pool_5aa, n)
                    for seq, set_num, attempt in triplication(
                        raw, patches, attempts, remove_prob
                    ):
                        yield {
                            'source':        'recombination',
                            'motif_type':    '5aa',
                            'motif_pattern': raw,
                            'n_repeats':     n,
                            'target_length': target,
                            'set':           set_num,
                            'attempt':       attempt,
                            'sequence':      seq,
                            'actual_length': len(seq),
                        }

        # ── 10AA recombination ─────────────────────────────────
        if pool_10aa:
            for n in valid_repeat_counts(10, target, tol):
                for _ in range(n_samples):
                    raw = _sample_combination(pool_10aa, n)
                    for seq, set_num, attempt in triplication(
                        raw, patches, attempts, remove_prob
                    ):
                        yield {
                            'source':        'recombination',
                            'motif_type':    '10aa',
                            'motif_pattern': raw,
                            'n_repeats':     n,
                            'target_length': target,
                            'set':           set_num,
                            'attempt':       attempt,
                            'sequence':      seq,
                            'actual_length': len(seq),
                        }

        # ── Mixed 5AA + 10AA recombination ────────────────────
        if pool_5aa and pool_10aa:
            for (a, b) in mixed_configs.get(target, []):
                for _ in range(n_samples):
                    motif_list = (
                        [random.choice(pool_5aa)  for _ in range(a)] +
                        [random.choice(pool_10aa) for _ in range(b)]
                    )
                    random.shuffle(motif_list)
                    raw = ''.join(motif_list)
                    for seq, set_num, attempt in triplication(
                        raw, patches, attempts, remove_prob
                    ):
                        yield {
                            'source':        'recombination',
                            'motif_type':    f'mixed_{a}x5aa_{b}x10aa',
                            'motif_pattern': raw,
                            'n_repeats':     a + b,
                            'target_length': target,
                            'set':           set_num,
                            'attempt':       attempt,
                            'sequence':      seq,
                            'actual_length': len(seq),
                        }


# ================================================================
#  OUTPUT & DEDUPLICATION
# ================================================================

FIELDNAMES = [
    'source',
    'motif_type',
    'motif_pattern',
    'n_repeats',
    'target_length',
    'set',
    'attempt',
    'sequence',
    'actual_length',
]


class DeduplicatorHashBased:
    """
    Memory-efficient deduplication using 128-bit MD5 hashes.
    At 20 M sequences → ~320 MB RAM.  Full-string storage → ~3–10 GB.
    Hash collisions are astronomically unlikely at this scale.
    """
    def __init__(self):
        self._seen: set = set()

    def is_new(self, seq: str) -> bool:
        h = hashlib.md5(seq.encode()).digest()
        if h in self._seen:
            return False
        self._seen.add(h)
        return True

    def __len__(self):
        return len(self._seen)


def _stream_and_write(
    generator: Generator[Dict, None, None],
    writer:    csv.DictWriter,
    dedup:     DeduplicatorHashBased,
    batch_size: int,
    label:     str,
) -> Tuple[int, int]:
    """Stream rows from generator → deduplicate → write to CSV in batches."""
    written = dupes = 0
    batch: List[Dict] = []
    t0 = time.time()

    for row in generator:
        if dedup.is_new(row['sequence']):
            batch.append(row)
            written += 1
        else:
            dupes += 1

        if len(batch) >= batch_size:
            writer.writerows(batch)
            batch.clear()
            elapsed = time.time() - t0
            print(
                f"\r  {label}: {written:>10,} written  "
                f"{dupes:>8,} dupes  "
                f"({elapsed:.0f}s)   ",
                end='', flush=True,
            )

    if batch:
        writer.writerows(batch)

    elapsed = time.time() - t0
    print(
        f"\r  {label}: {written:>10,} written  "
        f"{dupes:>8,} dupes  "
        f"({elapsed:.1f}s)   "
    )
    return written, dupes


# ================================================================
#  MAIN PIPELINE
# ================================================================

def run() -> None:
    if CONFIG['random_seed'] is not None:
        random.seed(CONFIG['random_seed'])

    _init_pools()

    tgt      = CONFIG['target_lengths']
    tol      = CONFIG['length_tolerance']
    sc       = CONFIG['sample_counts']
    attempts = CONFIG['set_attempts']
    rp       = CONFIG['set3_remove_prob']
    outfile  = CONFIG['output_file']
    bsz      = CONFIG['write_batch_size']

    sep = '=' * 62
    print(f"\n{sep}")
    print("  IDP Solubility Tag Sequence Generator")
    print(sep)

    # ── Stage 1: Motif pools ──────────────────────────────────
    print("\n[Stage 1]  Motif Generation")

    print("  Generating KPEV+ …")
    pool_kpev  = gen_kpev()
    print(f"            → {len(pool_kpev)} sequence")

    print("  Generating Motif 2 (4 AA) …")
    pool_4aa   = gen_motif2(sc['motif2'])
    print(f"            → {len(pool_4aa):,} sequences")

    print("  Generating 5 AA motifs …")
    pool_5aa   = gen_5aa_motifs(sc['5aa'])
    print(f"            → {len(pool_5aa):,} sequences")

    print("  Generating 10 AA motifs  (largest step) …")
    pool_10aa  = gen_10aa_motifs(sc['10aa'])
    print(f"            → {len(pool_10aa):,} sequences")

    patches    = gen_hydrophobic_patches()
    print(f"  Hydrophobic patches      → {len(patches)} patches")

    print(f"\n  Motif pool summary:")
    print(f"    KPEV+   : {len(pool_kpev):>10,}")
    print(f"    Motif 2 : {len(pool_4aa):>10,}")
    print(f"    5 AA    : {len(pool_5aa):>10,}")
    print(f"    10 AA   : {len(pool_10aa):>10,}")
    print(f"    Patches : {len(patches):>10,}")

    # ── Stage 2 & 3: Sequence generation ──────────────────────
    print(f"\n[Stage 2]  Repeat Algorithm")
    print(f"[Stage 3]  Recombination Algorithm")
    print(f"\n  Output → {outfile}\n")

    dedup        = DeduplicatorHashBased()
    total_written = 0
    total_dupes   = 0

    t_global = time.time()

    with open(outfile, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()

        # Repeat — KPEV+
        w, d = _stream_and_write(
            repeat_algorithm(pool_kpev, 'kpev', tgt, tol, patches, attempts, rp),
            writer, dedup, bsz, "KPEV+ repeats    ",
        )
        total_written += w; total_dupes += d

        # Repeat — Motif 2
        w, d = _stream_and_write(
            repeat_algorithm(pool_4aa, 'motif2', tgt, tol, patches, attempts, rp),
            writer, dedup, bsz, "Motif2 repeats   ",
        )
        total_written += w; total_dupes += d

        # Repeat — 5 AA
        w, d = _stream_and_write(
            repeat_algorithm(pool_5aa, '5aa', tgt, tol, patches, attempts, rp),
            writer, dedup, bsz, "5AA repeats      ",
        )
        total_written += w; total_dupes += d

        # Repeat — 10 AA  (largest batch)
        w, d = _stream_and_write(
            repeat_algorithm(pool_10aa, '10aa', tgt, tol, patches, attempts, rp),
            writer, dedup, bsz, "10AA repeats     ",
        )
        total_written += w; total_dupes += d

        # Recombination
        w, d = _stream_and_write(
            recombination_algorithm(
                pool_4aa, pool_5aa, pool_10aa,
                tgt, tol, patches, attempts, rp,
                sc['recombination'],
            ),
            writer, dedup, bsz, "Recombination    ",
        )
        total_written += w; total_dupes += d

    elapsed = time.time() - t_global
    fsize   = os.path.getsize(outfile)

    print(f"\n{sep}")
    print("  COMPLETE")
    print(sep)
    print(f"  Unique sequences written : {total_written:>12,}")
    print(f"  Duplicates removed       : {total_dupes:>12,}")
    print(f"  Total runtime            : {elapsed:>11.1f} s")
    print(f"  Output file              : {outfile}")
    print(f"  File size                : {fsize/1e6:>11.1f} MB")
    print(sep + "\n")


if __name__ == '__main__':
    run()
