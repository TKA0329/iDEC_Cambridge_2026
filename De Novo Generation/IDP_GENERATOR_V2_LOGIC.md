8888# IDP Generator v2 — Logic Documentation

## Overview
**Purpose:** In silico de novo generation of Intrinsically Disordered Protein (IDP) solubility tag candidate sequences.

**Output:** CSV file containing generated sequences with full provenance metadata (source, motif pattern, set type, attempt number, target length, actual length).

**Key Constraint:** Sequences are limited to ≤60 amino acids (sequences of exactly 60 AA are included; only sequences of 61 AA or longer are excluded) to remain experimentally tractable.

---

## 1. Original Design Prompt

The generator answers a fundamental biological question:

> **"How can we computationally generate novel IDP solubility tag sequences that maintain biological properties (charge patterns, hydrophobic-hydrophilic balance) while achieving target lengths?"**

Key assumptions:
- **Charge balance matters** → positive (K, R) and negative (E, D) residues are weighted
- **Hydrophobic content** → explicitly controlled via patches (G, A, V, I, L)
- **Polar residues** → solubility anchors (S, T, Q, N, H)
- **Length is parametric** → target ranges (25, 50 aa) with ±3 aa tolerance
- **Modularity enables scale** → building blocks (motifs) repeated/recombined to reach targets

---

## 2. Probing Questions & Design Decisions

The generator addresses these questions through its architecture:

### Q1: How do we encode known biological patterns?
**A:** Use fixed anchor motifs (KPEV+) and semi-random motifs (Motif 2, 5AA, 10AA) that incorporate weighted amino acid pools.

### Q2: How do we generate diversity within the charge/hydrophobic constraints?
**A:** Three generation strategies:
1. **Repeat Algorithm** — repeat a motif n times (simple, predictable)
2. **Recombination Algorithm** — mix different motif types (higher diversity)
3. **Set Triplication** — apply post-processing transformations to '+' placeholders

### Q3: What role does the '+' (plus) symbol play?
**A:** It's a **post-generation placeholder** for hydrophobic patches. The triplication step transforms it:
- **Set 1 (deterministic):** Remove all '+' → purely charge-polar backbone
- **Set 2 (random):** Replace each '+' with a randomly chosen hydrophobic patch → high hydrophobicity
- **Set 3 (probabilistic):** Each '+' is removed (75% default) or replaced by a patch (25%)
  (3 independent attempts)

### Q4: Why use weighted sampling for amino acids?
**A:** To bias generation toward biologically relevant sequences:
- Positive/negative: weight 1:1 (balanced charge)
- Non-polar set 1: P weighted 3× over G (Pro-rich motifs)
- Non-polar set 2: G weighted 2×, allowing V/I/L balance
- Polar: S/T/Q weighted 2×, N weighted 1×, H weighted 2× (conditional charge)

### Q5: How do we avoid combinatorial explosion?
**A:** Four strategies:
1. **Motif pools** — pre-generate a limited set of motifs (10–405 per type)
2. **Sampling** — take n_samples per configuration (10 for 4AA/5AA/recombination, 3 for 10AA)
3. **Mixed slot configs** — not all (a, b) pairs; only those that fit length tolerance
4. **Deduplication** — MD5 hash-based (memory-efficient, astronomically low collision risk)

### Q6: Why separate repeat and recombination?
**A:**
- **Repeat**: Simple, fast, predictable; good for baseline sequences
- **Recombination**: Higher complexity; explores novel combinations across motif types
- Both feed into the same triplication logic for consistency

### Q7: Why is hard cap filtering necessary?
**A:** To ensure all sequences remain experimentally feasible. Sequences >60 aa are excluded regardless of other parameters. Sequences of exactly 60 AA are retained.

---

## 3. Core Logic: Three-Stage Pipeline

### Stage 1: Motif Generation

Generates the building blocks from which all sequences are constructed.

#### 1.1 — Motif 1 (KPEV+)
- **Type:** Fixed biological anchor from Smt3 solubility tag
- **Unique sequences:** 1 (no sampling)
- **Structure:** K–P–E–V–+
- **Role:** Biological reference baseline
- **Note:** KPEV+ is passed only to the repeat algorithm and is **excluded from all recombination categories by design**.

#### 1.2 — Motif 2 (4 AA motifs)
- **Pattern:** (pos)–P–(neg)–np2–+
- **Sampling:** Up to 10 unique sequences (theoretical maximum is 20; if a higher sample count is requested, a warning is printed and the achievable unique count is returned instead)
- **Logic:**
  - Position 1: K or R
  - Position 2: Fixed P (Pro-rich)
  - Position 3: E or D
  - Position 4: G, A, V, I, or L
  - Position 5: Fixed + (placeholder)
- **Result:** Up to 10 sequences with controlled AA distribution

#### 1.3 — 5 AA motifs
- **Structure:** Each motif is **6 characters total** — 5 amino acid positions followed by a terminal `+` placeholder.
- **Pattern:** Exactly 1 non-polar (np1, position varies across positions 0–4) + 4 slots from {neg, pos, polar}, with `+` always appended as a 6th terminal character
- **Combinatorial:**
  - np1 position: 5 options
  - 4 remaining slots: 3^4 = 81 category combinations
  - Total patterns: 5 × 81 = 405
  - With n_samples=10 per pattern: ~4,050 unique sequences (theoretical upper bound; actual unique count will be somewhat lower due to random sampling collisions within the set-based deduplication)
- **Logic:** Each slot samples independently from its category

#### 1.4 — 10 AA motifs
- **Structure:** Each motif is **11 characters total** — 10 amino acid positions followed by a terminal `+` placeholder.
- **Pattern:** Exactly 1 non-polar (np1, position varies across positions 0–9) + 9 slots from {neg, pos, polar}, with `+` always appended as an 11th terminal character
- **Combinatorial:**
  - np1 position: 10 options
  - 9 remaining slots: 3^9 = 19,683 category combinations
  - Total patterns: 10 × 19,683 = 196,830
  - With n_samples=3 per pattern (reduced to manage scale): ~590,490 unique sequences (theoretical upper bound; actual unique count will be lower due to sampling collisions)
- **Progress tracking:** Real-time ETA display every 20K patterns

#### 1.5 — Hydrophobic patches
- **Definition:** All 25 ordered 2-AA combinations of {G, A, V, I, L}
- **Use:** Replacements for '+' symbols post-generation
- **Unique sequences:** 5 × 5 = 25 (GG, GA, GV, … LL)

### Stage 2: Repeat Algorithm

Expands each motif by repeating it n times to reach target lengths.

**Pseudocode:**
```
for each motif in pool:
    for each target_length in targets:
        for each valid_repeat_count n:
            raw_sequence = motif × n
            for each (set_num, resolved_seq) in triplication(raw, patches):
                if length(resolved_seq) ≤ cap:
                    yield {source: 'repeat', ...}
```

**Key steps:**
1. Calculate `aa_len` = motif length − 1 (excludes the single terminal `+`). This assumes each individual motif has exactly one `+`, which holds for all current motif types.
2. Find all n where: `(target - tolerance) ≤ n × aa_len ≤ (target + tolerance)`
3. Repeat motif n times: `KPEV+KPEV+KPEV+`
4. Apply triplication (see below)
5. Filter results under hard cap

**Example:**
- Motif: KPEV+ (5 chars total; 4 AA core excluding +)
- aa_len = 4
- Target: 25 aa ± 3 → [22, 28]
- Valid repeats: n ∈ {6, 7} → sequences of length 24 or 28
- CSV row example: `repeat, kpev, KPEV+, 7, 25, 1, 1, KPEVKPEVKPEVKPEVKPEVKPEVKPEV, 28`

### Stage 3: Recombination Algorithm

Mixes different motif types and samples random combinations to reach targets.

**Four recombination categories:**

1. **4AA-only recombination**
   - Recombine Motif 2 sequences with themselves
   - For each valid repeat count n (based on 4 AA core):
     - Sample n motifs randomly from pool_4aa
     - Apply triplication

2. **5AA-only recombination**
   - Recombine 5AA motifs with themselves
   - For each valid repeat count n (based on 5 AA core):
     - Sample n motifs from pool_5aa
     - Apply triplication

3. **10AA-only recombination**
   - Recombine 10AA motifs with themselves
   - For each valid repeat count n (based on 10 AA core):
     - Sample n motifs from pool_10aa
     - Apply triplication

4. **Mixed 5AA + 10AA recombination**
   - **Configuration:** For each target length, find all (a, b) pairs where:
     - a > 0, b > 0 (both required — pure-5AA and pure-10AA are handled by their own categories above)
     - `(target - tol) ≤ a × 5 + b × 10 ≤ (target + tol)`
   - **Sampling:** For each (a, b) pair:
     - Sample a random 5AA motifs
     - Sample b random 10AA motifs
     - Shuffle the list → randomized interleaving
     - Apply triplication

**Note:** KPEV+ (Motif 1) is excluded from all recombination categories and is processed only by the repeat algorithm.

**Example:**
- Target: 50 aa ± 3 → [47, 53]
- Possible configs: (7, 2) → 35 + 20 = 55 ✓, (9, 1) → 45 + 10 = 55 ✓
- For each config: sample 10 random combinations per config

### Triplication: Set 1 / 2 / 3 Processing

Each raw sequence (from repeat or recombination) undergoes **deterministic and stochastic transformations** to expand the space:

**Set 1 (Deterministic, 1 output per input):**
```
Remove all '+' symbols
Example: KPEV+RPEG+ → KPEVREPEG
```

**Set 2 (Fully random, `attempts` outputs per input, default=3):**
```
Replace every '+' with a random hydrophobic patch (25 options)
Example iteration 1: KPEV+RPEG+ → KPEVGIREPEGVI
Example iteration 2: KPEV+RPEG+ → KPEVAAREPEGVV
```

**Set 3 (Probabilistic, `attempts` outputs per input, default=3):**
```
Each '+' is independently:
  • Removed with probability remove_prob (default 0.75) → 75% removal
  • Replaced with patch with probability (1 - remove_prob) (default 0.25) → 25% replacement
Example iteration 1: KPEV+RPEG+ → KPEVREPEG    (both removed)
Example iteration 2: KPEV+RPEG+ → KPEVGIREPEGAA (both replaced)
Example iteration 3: KPEV+RPEG+ → KPEVGIREPEG   (first replaced, second removed)
```

**Output:**
- Set 1: 1 unique sequence per raw
- Set 2: up to 3 sequences per raw (with random patch choices)
- Set 3: up to 3 sequences per raw (with probabilistic removal/replacement)
- **Total triplication multiplier:** 1 + attempts + attempts = 1 + 2 × 3 = 7× per raw sequence

---

## 4. Final Output Logic

### 4.1 — Deduplication

**Strategy:** MD5 hash-based deduplication (memory-efficient)

**Why MD5 hashes instead of full string storage?**
- 20M sequences: ~320 MB RAM (MD5) vs 3–10 GB (full strings)
- Hash collision probability at 20M items: ~1 in 10^30 (astronomically unlikely)
- Trade-off: Deterministic collision detection (no false positives with proteins)

**Process:**
```
for each generated sequence:
    hash = MD5(sequence.encode())
    if hash not in seen_hashes:
        add to output rows
        seen_hashes.add(hash)
    else:
        increment duplicate counter
```

**Real-time progress tracking:**
- Updates every 10,000 sequences
- Shows: sequences kept, duplicates removed, elapsed time, running rate

### 4.2 — Sorting

All unique sequences are sorted by `actual_length` (ascending):
```
all_rows.sort(key=lambda r: r['actual_length'])
```

**Rationale:**
- Length is the primary biological parameter
- Sorting enables efficient downstream filtering (e.g., by length ranges)
- Reproducible output order across runs

### 4.3 — CSV Writing

**Output file:** `generated_sequences.csv` (configurable)

**Columns (FIELDNAMES):**
```
source              — 'repeat' or 'recombination'
motif_type          — 'kpev', 'motif2', '5aa', '10aa', or 'mixed_{a}x5aa_{b}x10aa'
motif_pattern       — Raw sequence(s) before triplication (with '+' placeholders)
n_repeats           — Number of times motif(s) were repeated/combined.
                      For mixed configurations, n_repeats = a + b (total motif slots),
                      not a simple repeat count.
target_length       — Target length aimed for [25, 50]
set                 — Triplication set (1, 2, or 3)
attempt             — Attempt number within set (1 for Set 1; 1–N for Sets 2 & 3)
sequence            — Final resolved sequence (no '+' symbols)
actual_length       — len(sequence)
```

**Write process:**
1. Open file in write mode
2. Write header row (FIELDNAMES)
3. Write all rows (deduped, sorted) in a single pass
4. Close file

**Note:** All rows are held in memory prior to writing (required for sorting by length). The `write_batch_size` parameter has been removed — batch flushing is not implemented and offers no benefit at this stage since all rows are already in memory.

### 4.4 — Final Report

Terminal output after completion:
```
================================================================
  IDP Solubility Tag Sequence Generator
================================================================
  Target lengths : [25, 50]  (±3 aa)
  Hard cap       : 60 aa

[Stage 1]  Motif Generation
  …pool generation …

[Stage 2]  Repeat Algorithm
[Stage 3]  Recombination Algorithm

  Sorting X,XXX,XXX sequences by actual_length …
  Writing to generated_sequences.csv …

================================================================
  COMPLETE
================================================================
  Unique sequences written : X,XXX,XXX
  Duplicates removed       : XXX,XXX
  Total runtime            : XXX.X s
  Output file              : generated_sequences.csv
  File size                : XXX.X MB
================================================================
```

---

## 5. Configuration Parameters

All tunable parameters reside in `CONFIG` dict:

```python
CONFIG = {
    'sample_counts': {
        'motif2':        10,    # Motif 2 samples per pattern
        '5aa':           10,    # 5 AA samples per (position, combo)
        '10aa':           3,    # 10 AA samples per (position, combo) — reduced for scale
        'recombination': 10,    # Samples per slot-config × target × motif-category
    },
    'ratios': {
        'positive':  {'K': 1, 'R': 1},
        'negative':  {'E': 1, 'D': 1},
        'np1':       {'G': 1, 'P': 3},
        'np2':       {'G': 2, 'A': 1, 'V': 1, 'I': 1, 'L': 1},
        'polar':     {'S': 2, 'T': 2, 'Q': 2, 'N': 1, 'H': 2},
    },
    'target_lengths':   [25, 50],
    'length_tolerance': 3,
    'length_hard_cap':  60,
    'set_attempts':      3,        # Repeat count for Sets 2 & 3
    'set3_remove_prob':  0.75,     # Probability of removing '+' in Set 3
    'output_file':      'generated_sequences.csv',
    'random_seed': 42,             # None for different run each time
}
```

**Reproducibility note:** Setting `random_seed` guarantees identical output only when all CONFIG parameters are also identical. Changing any parameter (including sample counts) shifts the global random state and will produce different Stage 2/3 outputs even with the same seed. Reproducibility is therefore global across the full run, not per-stage.

**Tuning implications:**
- Increase `sample_counts` → more sequences, longer runtime
- Increase `target_lengths` → more configurations to explore
- Increase `set_attempts` → 7× triplication multiplier becomes (1 + 2×N)×
- Increase `length_hard_cap` → more sequences survive filtering (note: sequences of exactly the cap value are retained)
- Adjust ratios → change AA probabilities (e.g., more Pro, less Gly)

---

## 6. Data Flow Summary

```
┌─────────────────────────────────┐
│  Stage 1: Motif Pools           │
│  • KPEV+ (1) — repeat only      │
│  • Motif 2 (≤10)                │
│  • 5 AA (≤4,050)                │
│  • 10 AA (≤590,490)             │
│  • Patches (25)                 │
└──────────────┬──────────────────┘
               │
         ┌─────┴─────┐
         │           │
    ┌────▼────┐  ┌───▼──────────────────────┐
    │ Repeat  │  │ Recombination            │
    │ Algo    │  │ (4AA/5AA/10AA/mixed)     │
    │(all     │  │ (KPEV+ excluded)         │
    │ pools)  │  │                          │
    └────┬────┘  └───┬──────────────────────┘
         │           │
         └─────┬─────┘
               │
         ┌─────▼─────┐
         │Triplication│
         │(Set 1/2/3) │
         └─────┬─────┘
               │
         ┌─────▼──────────────┐
         │ Deduplication      │
         │ (MD5 hash-based)   │
         └─────┬──────────────┘
               │
         ┌─────▼──────────────┐
         │ Sort by Length     │
         └─────┬──────────────┘
               │
         ┌─────▼──────────────┐
         │ Write CSV          │
         │ (with metadata)    │
         └────────────────────┘
```

---

## 7. Example Execution Flow

**Input:** (None — fully procedural)

**Execution steps:**
1. Initialize random seed (42)
2. Build amino acid pools (weighted sampling)
3. Generate motif pools (Stage 1)
4. Repeat algorithm on KPEV+, Motif 2, 5AA, 10AA → collector
5. Recombination algorithm (4AA/5AA/10AA/mixed, KPEV+ excluded) → collector
6. Collect results with deduplication and progress tracking
7. Sort all rows by `actual_length`
8. Write to CSV

**Example output snippet (first 5 rows after sorting):**
```
source,motif_type,motif_pattern,n_repeats,target_length,set,attempt,sequence,actual_length
repeat,kpev,KPEV+,6,25,1,1,KPEVKPEVKPEVKPEVKPEVKPEV,24
repeat,kpev,KPEV+,7,25,2,1,KPEVKPEVKPEVKPEVKPEVKPEVKPEVGA,30
repeat,motif2,RPEG+KPEV+…,6,25,1,1,RPEGKPEV…,24
recombination,mixed_2x5aa_2x10aa,GPSQ+KPDE+GPSN+KPDA+…,4,25,3,2,GPSQKPDEGSQKPDAA,25
…
```

**Note on KPEV+ repeat counts and target lengths:** KPEV+ has a 4 AA core (excluding `+`). For target=25 ± 3 ([22, 28]), valid repeat counts are n=6 (24 AA) and n=7 (28 AA). These rows are recorded with `target_length=25`. For target=50 ± 3 ([47, 53]), valid repeat counts are n=12 (48 AA) and n=13 (52 AA).

---

## 8. Computational Complexity

| Stage | Operation | Time Complexity | Notes |
|-------|-----------|-----------------|-------|
| Stage 1 | Motif generation | O(n_samples × n_patterns) | 10 AA slowest (196K patterns) |
| Stage 2 | Repeat algorithm | O(n_motifs × n_targets × n_repeats × triplication_factor) | Predictable, fast |
| Stage 3 | Recombination | O(n_targets × n_configs × n_samples × triplication_factor) | Stochastic; largest contributor |
| Dedup | MD5 hashing | O(n_total_sequences) | Linear, memory-efficient |
| Sort | Quicksort | O(n_unique × log(n_unique)) | After dedup |
| Write | CSV I/O | O(n_unique) | Single-pass write after sort |

**Typical runtime:** 10–60 seconds (varies by configuration and hardware)

---

## 9. Key Design Decisions & Justifications

| Decision | Rationale |
|----------|-----------|
| Fixed KPEV+ anchor | Proven biological solubility tag; provides baseline |
| KPEV+ repeat-only | Excluded from recombination by design to preserve its fixed-anchor role |
| Weighted amino acid pools | Stochastic generation with deterministic bias |
| '+' placeholder strategy | Decouples motif design from post-processing; enables Set 1/2/3 diversity |
| Hard cap at 60 aa (inclusive) | Experimental tractability; sequences of exactly 60 AA are retained |
| Length tolerance ±3 | Allows flexibility without overfitting to exact target |
| MD5 deduplication | Memory-efficient; collision risk negligible at this scale |
| Sort by length | Primary biological parameter; enables downstream filtering |
| Triplication (3 sets) | Expands sequence diversity without explicit motif multiplication |
| Recombination interleaving | Random shuffle explores more sequence space than simple concatenation |
| Global random seed | Reproducibility requires identical CONFIG parameters; changing any parameter shifts the global random state for all downstream stages |

---

## 10. Future Enhancements

Possible extensions to the generator:

1. **Per-stage seeding** — seed each stage independently (e.g. `random.seed(base_seed + stage_number)`) to allow parameter tuning without shifting all downstream random state
2. **GPU-accelerated sampling** — parallel motif generation
3. **Constraint-based generation** — enforce specific charge ratios or hydrophobicity scores during generation rather than post-hoc filtering
4. **Predicted property pre-filtering** — calculate GRAVY, pI, secondary structure during generation
5. **Sequence similarity clustering** — remove near-duplicates (Levenshtein distance < threshold)
6. **Theoretical maximum enforcement in Motif 2** — compute the 20-sequence alphabet ceiling upfront and cap the requested sample count before entering the generation loop, rather than relying on the max_tries warning