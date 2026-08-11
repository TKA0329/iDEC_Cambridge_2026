# Degenerate-motif proportion sweep pipeline

`run_pipeline.py` is the new orchestrator. It sweeps the `--proportion`
parameter of `filter_sequences_pt2.py` and, for every value, runs the full
chain:

```
filter_sequences_pt2.py  -->  reverse_translation_V5_StopSafe.py  -->  sample_degenerate.py  -->  batch_analyzer.py  -->  filter_sequences.py
   (build motif)                 (--mode direct, degenerate DNA)         (sample + translate)      (ProtParam/CIDER)      (--non-interactive)
```

...and writes one `summary.csv` row per proportion (`proportion,
degenerate_sequence, encoded_space, samples_generated, rows_analyzed,
rows_before_filter, rows_after_filter, pct_pass`), plus a full folder of
every intermediate file per proportion, so you can inspect/re-use any step.

## Quick start

```bash
python run_pipeline.py your_filtered_set.fasta my_run \
    --prop-min 0.1 --prop-max 0.9 --prop-step 0.1 \
    --num-samples 100000 --seed 42 \
    --gravy-max 0.5 --instability-max 40 \
    --ncpr-min -0.3 --ncpr-max 0.3
```

All values shown above except the fasta/output-dir are already the defaults
(0.1-0.9 step 0.1, 100000 samples, seed 42) — the filter thresholds (`--gravy-max`
etc.) are the one thing you must actively decide, since there's no sane
universal default for "good" GRAVY/instability/etc. cutoffs. Run
`python run_pipeline.py --help` to see every flag and its default; every
default is overridable.

Output layout:

```
my_run/
  summary.csv                  <- one row per proportion, this is the headline result
  prop_0.1/
    variants.fasta              filter_sequences_pt2.py output (degenerate protein motif)
    pt2_report.tsv              filter_sequences_pt2.py --report
    degenerate.fasta            RT's resulting degenerate DNA sequence (NEW, see below)
    02_rt.log                   full RT terminal output (RT's own --log)
    samples.csv                 sample_degenerate.py output (protein sequences)
    analyzed.csv                batch_analyzer.py output (ProtParam/CIDER columns)
    filtered/
      filtered_analyzed.csv     final filtered rows + metadata comments
      length_summary_stats.csv
    0N_*.log                    orchestrator's own copy of each step's stdout
  prop_0.2/ ...
  ...
```

The best proportion *directly sampled* by the sweep is printed at the end
and is in `summary.csv` (sort by `pct_pass`). This is a brute-force grid
search, not a curve fit — per your note, that's intentional for v1. Once
you've seen the shape of the curve, narrow `--prop-min/--prop-max/--prop-step`
around the best point and re-run for a finer local search; a regression-based
optimizer (matching the style of RT's own `--bins-sweep`/`--stop-sweep`
tradeoff tables) can be added later if the grid search isn't precise enough.

## What was changed in the existing scripts, and why

1. **`sample_degenerate.py`** — this is the script you meant instead of
   `iupac_sampler.py` (which is no longer used in the pipeline; it's still
   here but not called by `run_pipeline.py`). Two edits:
   - `--seq-file` now tolerates FASTA input (strips `>` header lines), so it
     can read RT's `--out-seq` output directly. If the file has more than one
     FASTA record (i.e. RT was run with `--bins > 1`), only the first
     record's sequence is used, with a warning — this pipeline runs RT with
     `--bins 1` by default, so that's not a concern unless you deliberately
     raise `--rt-bins`.
   - The output CSV's sequence column header changed from `sequences` to
     `sequence` (singular) to match what `sequence_stats.py` (and therefore
     `filter_sequences.py`) actually requires. This was a pre-existing
     mismatch between the two scripts that would have made the final filter
     step fail with "Missing required columns: sequence".

2. **`reverse_translation_V5_StopSafe.py`** — added one new flag, `--out-seq
   FILE`. Previously the script only *printed* the resulting degenerate
   sequence(s); nothing was ever written to disk (only `--bins-sweep`/
   `--stop-sweep` populate `--csv`, and that's a comparison table across
   sweep steps, not the final answer). `--out-seq` writes the single-run
   result as FASTA (`>bin1`, `>bin2`, ... one record per bin) so it can be
   piped straight into `sample_degenerate.py`. It's a no-op (with a warning)
   during `--bins-sweep`/`--stop-sweep`, and only fires for the plain
   single-run path. Everything else in the script — grid/surface modes,
   the sweep modes, manual weights — is untouched, so you can still use them
   by hand exactly as before; the orchestrator just doesn't call them.

3. **`filter_sequences.py`** — added a `--non-interactive` mode with one CLI
   flag per threshold (`--gravy-max`, `--instability-max`, `--helix-sheet-max`,
   `--turn-fraction-min`, `--charge-sum-min`, `--charge-ph7-exclude-min`/`-max`,
   `--fcr-max`, `--ncpr-min`/`-max`, `--cider-kappa-max`), so the orchestrator
   can call it headlessly with the exact same thresholds at every proportion
   (apples-to-apples comparison). Run it with no flags and it behaves exactly
   as before — interactive prompts, unchanged. Also added a `pass_rate_pct`
   field to the metadata comments it writes, so downstream tools (the
   orchestrator, or you by hand) don't have to recompute it from
   rows_before/rows_after.

4. **`filter_sequences_pt2.py`, `batch_analyzer.py`, `protein_pipeline.py`,
   `sequence_stats.py`** — unchanged, used as-is (already fully CLI-driven).

## Notes / things worth knowing

- **Comparability across proportions**: the `pct_pass` numbers are only
  meaningful relative to each other because the filter thresholds are fixed
  once for the whole sweep. If you change any `--*-max`/`--*-min` flag
  between runs, don't compare `pct_pass` across those runs.
- **Sampling noise**: `sample_degenerate.py --unique` (the default) draws
  without replacement up to the degenerate sequence's total combinatorial
  space; if `encoded_space` for a proportion is smaller than `--num-samples`,
  you'll see a "could only generate N unique sequences" warning in that
  proportion's `03_sample.log`, and `samples_generated` in `summary.csv` will
  be less than requested. In that situation `pct_pass` is being computed over
  the *entire* enumerable space (a census, not a sample) so it's actually
  more reliable, not less.
- **Failure handling**: by default a failed proportion (e.g. RT can't find any
  stop-free codon covering enough sequences) doesn't kill the sweep — it's
  logged with `status = FAILED: ...` in `summary.csv` and the orchestrator
  moves on. Pass `--stop-on-error` if you'd rather it stop immediately.
- **Reproducibility**: `--seed` (default 42) is reused for RT's internal
  sampling/clustering and for `sample_degenerate.py`'s random draw. It's only
  passed to `filter_sequences_pt2.py` when `--random-base` is set (the default
  consensus-base mode is deterministic and doesn't use randomness).
