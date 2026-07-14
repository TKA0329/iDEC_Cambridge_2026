# cahs_analyzer_03.py 
* streaming scorer + filter for very large mutant libraries (built/tested for
10M-row, 1GB+ CSVs). Scores each sequence based on these properties: Pace & Scholtz helix propensity, Eisenberg hydrophobic moment, helical face occupancy, salt-bridge count, Pro/Gly hard filters and writes ONLY the rows that pass the thresholds -- nothing is held in memory across rows, so this scales to files far bigger than the RAM.

* Allows other biopython parameters to be switched on such as gravy, instability_index, aliphatic index, net charge at pH 7 and isoelectric point. 

KEY DIFFERENCE FROM cahs_analyzer_02:
- Old script: read all rows into a list, scored them, wrote all rows out.
  Fine for hundreds of rows; at 10M rows this can use many GB of RAM.
- This script: reads one row, scores it, checks filters, writes it (or
  doesn't), and immediately moves on. Memory use stays flat regardless of
  file size.

USAGE:
    python cahs_analyzer_03.py input.csv output.csv

INPUT CSV:
    Must have a column named "sequence" (case-insensitive). An optional
    "id" column is carried through.

OUTPUT CSV:
    Only rows that PASS all active filters, with all scored columns.
    A summary of how many rows were excluded by each criterion is printed
    to the terminal at the end (nothing extra is written to disk, so this
    won't blow up your output size).

# extract_glycine_motif.py
* Pulls out just the rows from a scored CSV where motif1_sequence contains at least one
Glycine, and adds columns describing WHERE in the motif that Gly sits.

* Adds:
    glycine_positions   -> 1-indexed position(s) of G within motif1_sequence,
                           comma-separated if more than one (e.g. "5" or "5,12")
    glycine_region       -> front / middle / back, based on which third of
                           motif1's length each G falls into. If a sequence
                           has multiple Gly in different regions, all
                           regions are listed (e.g. "front,middle")

USAGE:
    python extract_glycine_motifs.py scored_input.csv glycine_only_output.csv


# per_pos_frequency.py
* Computes the numeric version of your WebLogo: for each position in
motif1_sequence across your whole shortlisted CSV, counts how often each
amino acid appears, and reports the smallest set of residues needed to
cover a given fraction (default 90%) of survivors at that position.

* This is meant to run on your ROUND 2 shortlist (the CSV with a
motif1_sequence column, e.g. output of motif1_filter_pipeline.py), since
that's the data that should drive your actual degenerate-codon choices per
hotspot position.

USAGE:
    python position_frequency_analyzer.py shortlist.csv position_frequencies.csv

OUTPUT CSV columns:
    position, amino_acid, count, percentage, cumulative_percentage

Also prints, per position, a short summary: how many distinct residues
appear, and the minimal residue set needed to cover >= COVERAGE_THRESHOLD
of all survivors at that position (useful for picking a degenerate codon
narrower than full NNS where the data supports it).

