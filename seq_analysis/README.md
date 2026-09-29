# Comparing two variant compositions

This script follows four steps: read and align counts; calculate Pearson's chi-square; simulate null contingency tables with fixed margins; estimate an upper-tail p-value. It uses exact amino-acid strings as categories, with no additional error correction or clustering. The uploaded files contain 33 and 19 raw distinct variants, matching the raw counts rather than the clustered summaries in the earlier PDF.

## Run

Requires Python 3, NumPy and SciPy >=1.10. Tested with NumPy 2.3.5 and SciPy 1.17.0.

```bash
python -m pip install numpy "scipy>=1.10"
python compare_compositions.py "ssb_aa_variants 0.csv" "ssb_aa_variants 40k.csv" --simulations 100000 --seed 42 --outdir results
```

Both CSVs must have columns `aa_variant,count`. Options `--variant-col` and `--count-col` support other headers. Repeated identical sequence rows are summed. Counts must be nonnegative integers. Missing sequences in the other full export receive count zero; a top-ten-only export is not suitable. Columns with zero counts in both samples are discarded. No pseudocounts are added and totals need not be equal.

## What each step means

1. **Inputs and observed statistic.** Align the union of sequence identities into a 2 x K table. Expected counts are row_total * column_total / grand_total. Compute X² = sum((observed - expected)² / expected). No Yates correction is applied.
2. **Random tables.** The null hypothesis is equal variant proportions in the two sampled populations. Conditioning on the observed totals, randomise allocation of pooled variant counts between the two samples. SciPy's `random_table.rvs` with Patefield's algorithm samples the appropriate fixed-margin null distribution. It does not select all possible tables with equal probability. Both each sample's total reads and each variant's pooled count stay fixed.
3. **Simulated statistics.** Calculate the same Pearson statistic for each table. Expected counts stay the same because margins stay fixed. Simulations run in batches to limit memory use.
4. **P-value.** Count b simulated statistics at least as large as observed, including ties. With B simulations, report (b+1)/(B+1). This Monte Carlo correction avoids a zero p-value. The smallest reportable value is 1/(B+1); reaching it does not establish the actual tail probability at that precision. The JSON includes an exact binomial interval for Monte Carlo tail-probability uncertainty, not biological uncertainty.

This is a conditional Monte Carlo test using Pearson's statistic. It uses neither the asymptotic chi-square reference distribution nor Fisher's probability-based ordering of tables.

## Outputs

- `results/summary.json`: observed statistic, totals, expected-count diagnostics, simulated exceedances, corrected Monte Carlo p-value, seed, package versions, input SHA-256 hashes, and descriptive Cramér's V.
- `results/aligned_counts.csv`: complete aligned variant counts, frequencies, expected counts and each variant's contribution to X². Contributions are descriptive, not individual variant p-values.

## Result for the uploaded SSB files

- Untreated: 436 accepted reads, 33 exact amino-acid variants.
- Selected: 291 accepted reads, 19 exact amino-acid variants.
- Union: 51 variants; only one exact variant is shared.
- Pearson X²: 682.3544239360033.
- Expected cells below five: 87/102 (85.3%).
- Simulations: 100,000; seed: 42.
- Simulated statistics >= observed: 0.
- Corrected Monte Carlo p-value: 1/100001 = 0.000009999900001 (simulation floor).
- 95% binomial interval for the simulated tail probability: 0 to 0.0000368881. This quantifies simulation uncertainty only.
- Descriptive Cramér's V: 0.9688.

Suggested reporting: “The raw amino-acid variant compositions differed between the two sequenced samples (Pearson X² = 682.35; fixed-margin Monte Carlo test, 100,000 simulations; no simulated statistic equalled or exceeded the observed statistic; corrected Monte Carlo p = 1.00 x 10^-5, at the simulation resolution limit).”

## Interpretation and limitations

The result is conditional on the supplied accepted-read subset. It is not a test of repeatable selection effects across biological replicates and cannot resolve growth, preparation, amplification, QC or sequencing biases. Read independence and representative sampling remain assumptions. A large difference can be driven by different source populations or processing as well as selection.

These exports are unclustered: sequence errors and genuine rare variants are both categories. Repeat the analysis after final QC and, if used, a consistent cross-sample clustering procedure. Do not equate a good motif call with a viable biological construct. Use Batch 3's matching untreated sample rather than pooling unrelated untreated batches. Confirm sample lineage independently.

For several planned library comparisons, adjust the resulting global p-values for multiplicity. This script performs one global comparison and does not implement individual-variant tests.

## Verification

The implementation was checked against SciPy's Pearson statistic, fixed simulated row/column totals, seed reproducibility, and a fully enumerated small 2 x 2 null distribution. Input data and repository scripts were not modified.

Reference: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.random_table.html
