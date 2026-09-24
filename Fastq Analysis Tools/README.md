# Generic degenerate-sequence FASTQ analysis

This directory contains a reference-guided pipeline for measuring variants in sequencing reads from any DNA sequence library. It is not tied to a particular gene, tag, organism, or library name.

The analyzer compares:

- an original, non-degenerate DNA reference;
- a design sequence containing IUPAC ambiguity codes at permitted variable positions; and
- a FASTQ file containing sequencing reads.

The two FASTA records must be the same length and positionally aligned. The design sequence may contain `A`, `C`, `G`, `T` and standard IUPAC codes such as `R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V` and `N`.

## How variants are counted

A variant is the haplotype formed by the bases at IUPAC-degenerate positions only. Errors at fixed positions therefore do not create false unique variants. For each primary mapped read, the pipeline:

1. maps it to an unambiguous version of the design sequence;
2. calls each designed position using the read base and its quality;
3. rejects incomplete, low-quality, or design-invalid calls;
4. collapses identical valid haplotypes;
5. reconstructs the complete analyzed sequence for each haplotype; and
6. compares each reconstructed sequence with the original reference.

Protein outputs are optional. Use `--translate` only when the analyzed region is a coding sequence whose length is divisible by three.

## Installation

On Ubuntu:

```bash
sudo apt update
sudo apt install minimap2 samtools python3
```

No Python packages are required by `analyze_variants.py`.

## Input FASTA

The simplest input has two records. Record order is used by default:

```fasta
>original
ACGTACGTACGT
>design
ACGTACNTACGT
```

You can use any record names and select them explicitly with `--original-name` and `--degenerate-name`.

## Run the pipeline

From this directory:

```bash
python3 analyze_variants.py \
  --fasta reference_and_design.fasta \
  --fastq reads.fastq \
  --outdir analysis
```

The convenience wrapper has the same generic behavior:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh reference_and_design.fasta reads.fastq analysis
```

With named FASTA records:

```bash
./run_pipeline.sh reference_and_design.fasta reads.fastq analysis original design
```

To analyze only a region, provide 1-based inclusive coordinates. Coordinates refer to the full FASTA sequences:

```bash
python3 analyze_variants.py \
  --fasta reference_and_design.fasta \
  --fastq reads.fastq \
  --outdir region_analysis \
  --region-start 101 \
  --region-end 400
```

For a coding region, add `--translate`:

```bash
python3 analyze_variants.py \
  --fasta reference_and_design.fasta \
  --fastq reads.fastq \
  --outdir coding_analysis \
  --translate
```

Useful quality options are `--min-mapq` (default `20`), `--min-baseq` (default `8`) and `--min-fixed-baseq` (default `12`). Use `--preset` to select the appropriate minimap2 preset for the sequencing technology, for example `map-ont` or `sr`.

## Outputs

- `qc_summary.txt`: read counts, filtering counts, observed haplotypes, and design size.
- `variants.csv`: one row per valid nucleotide haplotype, including read count, frequency, reconstructed sequence, and changes from the original.
- `position_frequencies.csv`: base counts and frequencies at every designed position.
- `fixed_position_mismatches.csv`: QC for unexpected changes at fixed positions.
- `invalid_haplotypes.csv`: calls containing bases forbidden by the design.
- `read_assignments.tsv`: auditable status for every primary mapped read.
- `protein_variants.csv`: generated when `--translate` is used; synonymous DNA haplotypes are collapsed by encoded protein.
- `mapping_reference.fasta`, `aligned.sorted.bam`, and its index: mapping files.

`frequency` in the variant files is the number of reads assigned to that variant divided by all confidently assigned reads.

The `fastq_to_csv.py` utility is independent of the variant analyzer. It reports per-read length, base composition, quality, estimated error rate, and longest homopolymer:

```bash
python3 fastq_to_csv.py reads.fastq read_statistics.csv
```

## Interpretation

Read frequencies estimate the sequencing-library population, not necessarily the exact biological cell frequency. PCR, extraction, library preparation, and sequencing can introduce abundance bias. Comparisons between samples are most defensible when they use the same preparation and sequencing workflow.
