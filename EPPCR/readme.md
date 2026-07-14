# Version 2
Simulate an error-prone PCR mutant library.

* Error-prone PCR typically introduces largely random, independent point
substitutions at each nucleotide position. This simulation treats the
provided error rate as the PER-CYCLE per-base mutation probability, and
compounds it across --cycles to get the final per-base probability in the
amplified product:

    final_per_base_probability = 1 - (1 - per_cycle_rate) ** cycles

* For example, a per-cycle rate of 0.6% (0.006) over 25 cycles gives a final
per-base probability of ~14.1%, i.e. roughly 141 mutations per kb -- NOT
6 mutations per kb. If you already have a "final" rate in hand (e.g. from a
kit datasheet reporting mutations/kb after N cycles) either set --cycles 1
so no compounding occurs, or back-calculate the per-cycle rate first.

* The final mutation count for a molecule is drawn from a binomial distribution
over the sequence length, treating every position as an independent site that
may acquire a single substitution (positions are drawn without replacement,
so no site is hit twice in a single output molecule).

* Optionally, molecule-to-molecule rate variability can be introduced with a
gamma distribution, and/or a custom substitution spectrum (transition bias
per reference base) and per-position hotspot weighting can be supplied to
better approximate a specific polymerase's known error bias.

Example:
    python3 eppcr.py \
        --sequence ATGGCTGCTGCTTAA \
        --count 100000 \
        --error-rate 0.6 --units percent \
        --cycles 25 \
        --output mutants.csv \
        --summary mutations.tsv \
        --seed 42

# Version 1 
* Simulate an error-prone PCR mutant library.

* Error-prone PCR typically introduces largely random, independent point
substitutions at each nucleotide position. This simulation treats the provided
error rate as the final per-base mutation probability in the amplified product.

* For example, 0.6% means 0.006 mutations per nucleotide, or 6 mutations per kb.

* The final mutation count for a molecule is drawn from a binomial distribution
over the sequence length, treating every position as an independent site that
may acquire a single substitution.

* Optionally, molecule-to-molecule rate variability can be introduced with a
gamma distribution.

Example:
    python3 eppcr.py \
        --sequence ATGGCTGCTGCTTAA \
        --count 100000 \
        --error-rate 0.006 \
        --cycles 25 \
        --output mutants.csv \
        --summary mutations.tsv \
        --seed 42