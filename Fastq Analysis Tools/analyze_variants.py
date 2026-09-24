#!/usr/bin/env python3
"""
Reference-guided analysis of a degenerate sequence library.

Inputs
------
1. FASTA with two records:
    - original reference sequence
    - degenerate design sequence (IUPAC ambiguity codes)
2. FASTQ containing long-read sequencing results.

External tools
--------------
- minimap2
- samtools

Python dependencies
-------------------
Standard library only.

Main outputs
------------
qc_summary.txt
variants.csv
protein_variants.csv
position_frequencies.csv
read_assignments.tsv
invalid_haplotypes.csv
fixed_position_mismatches.csv
mapping_reference.fasta
aligned.sorted.bam (+ index)
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Call designed library haplotypes from long-read FASTQ data."
    )
    p.add_argument("--fasta", required=True,
                   help="Two-record FASTA containing original and degenerate references.")
    p.add_argument("--fastq", required=True, help="Sequencing FASTQ.")
    p.add_argument("--outdir", default="variant_analysis",
                   help="Output directory (default: variant_analysis).")
    p.add_argument("--original-name", default=None,
                   help="FASTA record name for the original sequence (default: first record).")
    p.add_argument("--degenerate-name", default=None,
                   help="FASTA record name for the degenerate design (default: second record).")
    p.add_argument("--region-start", type=int, default=None,
                   help="1-based inclusive sequence region start (default: 1).")
    p.add_argument("--region-end", type=int, default=None,
                   help="1-based inclusive sequence region end (default: end of sequence).")
    p.add_argument("--translate", action="store_true",
                   help="Translate the analyzed region and write protein variant outputs. "
                        "The region length must be divisible by three.")
    p.add_argument("--min-mapq", type=int, default=20,
                   help="Minimum mapping quality (default: 20).")
    p.add_argument("--min-baseq", type=int, default=8,
                   help="Minimum base quality at every designed variable position (default: 8).")
    p.add_argument("--min-fixed-baseq", type=int, default=12,
                   help="Minimum base quality for fixed-position mismatch QC (default: 12).")
    p.add_argument("--preset", default="map-ont",
                   help="minimap2 preset (default: map-ont).")
    p.add_argument("--skip-mapping", action="store_true",
                   help="Reuse aligned.sorted.bam already present in outdir.")
    return p.parse_args()


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence encountered before first header.")
                chunks.append(line)
        if name is not None:
            records[name] = "".join(chunks).upper()
    return records


def translate(seq):
    if len(seq) % 3:
        raise ValueError(f"Translated region length {len(seq)} is not divisible by 3.")
    return "".join(CODON_TABLE.get(seq[i:i+3], "X")
                   for i in range(0, len(seq), 3))


def resolve_for_mapping(degen_seq, original_seq, region):
    """
    Build an A/C/G/T-only mapping reference.
    At ambiguous analyzed positions, prefer the original base if it is allowed;
    otherwise choose a deterministic allowed base. Outside the analyzed region,
    replace any ambiguity with a deterministic allowed base.
    """
    chars = list(degen_seq)
    start, end = region

    for i, ch in enumerate(chars):
        if ch in "ACGT":
            continue
        if ch not in IUPAC:
            raise ValueError(f"Unsupported IUPAC symbol {ch!r} at position {i+1}.")
        allowed = IUPAC[ch]
        replacement = None
        if start <= i < end:
            ob = original_seq[i]
            if ob in allowed:
                replacement = ob
        if replacement is None:
            replacement = sorted(allowed)[0]
        chars[i] = replacement
    return "".join(chars)


def check_tools():
    missing = [x for x in ("minimap2", "samtools") if shutil.which(x) is None]
    if missing:
        raise SystemExit(
            "Missing required command-line tools: "
            + ", ".join(missing)
            + "\nOn Ubuntu install with:\n"
              "  sudo apt update\n"
              "  sudo apt install minimap2 samtools python3"
        )


def run_mapping(map_ref, fastq, bam, preset):
    sam = bam.with_suffix(".sam")
    print("[1/3] Mapping reads with minimap2...", file=sys.stderr)
    with open(sam, "w") as out:
        subprocess.run(
            ["minimap2", "-ax", preset, "--secondary=no", str(map_ref), str(fastq)],
            stdout=out, check=True
        )

    print("[2/3] Sorting/indexing BAM with samtools...", file=sys.stderr)
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    sam.unlink(missing_ok=True)


_cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_pairs(pos1, cigar, seq, qual):
    """
    Return dict: 0-based reference coordinate -> (base_or_dash, phred)
    for one SAM alignment.

    SAM POS is 1-based. SEQ/QUAL are already represented in alignment orientation
    as written in SAM/BAM, so CIGAR can be applied directly.
    """
    ref = pos1 - 1
    qry = 0
    obs = {}

    if qual == "*":
        qvals = [99] * len(seq)
    else:
        qvals = [ord(c) - 33 for c in qual]

    for n_str, op in _cigar_re.findall(cigar):
        n = int(n_str)

        if op in ("M", "=", "X"):
            for j in range(n):
                if qry + j < len(seq):
                    obs[ref + j] = (seq[qry + j].upper(), qvals[qry + j])
            ref += n
            qry += n
        elif op == "I":
            qry += n
        elif op in ("D", "N"):
            for j in range(n):
                obs[ref + j] = ("-", 0)
            ref += n
        elif op == "S":
            qry += n
        elif op in ("H", "P"):
            pass
        else:
            raise ValueError(f"Unexpected CIGAR operation: {op}")
    return obs


def sam_alignments_from_bam(bam):
    proc = subprocess.Popen(
        ["samtools", "view", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            yield {
                "qname": f[0],
                "flag": int(f[1]),
                "rname": f[2],
                "pos": int(f[3]),
                "mapq": int(f[4]),
                "cigar": f[5],
                "seq": f[9],
                "qual": f[10],
            }
    finally:
        proc.stdout.close()
        rc = proc.wait()
        if rc:
            raise RuntimeError(f"samtools view exited with code {rc}")


def nt_changes(original, variant):
    changes = []
    for i, (a, b) in enumerate(zip(original, variant), start=1):
        if a != b:
            changes.append(f"{a}{i}{b}")
    return changes


def aa_changes(original, variant):
    changes = []
    for i, (a, b) in enumerate(zip(original, variant), start=1):
        if a != b:
            changes.append(f"{a}{i}{b}")
    return changes


def main():
    args = parse_args()
    check_tools()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    records = read_fasta(args.fasta)
    if args.original_name not in records:
        raise SystemExit(
            f"Original record {args.original_name!r} not found. "
            f"Available: {', '.join(records)}"
        )
    if args.degenerate_name not in records:
        raise SystemExit(
            f"Degenerate record {args.degenerate_name!r} not found. "
            f"Available: {', '.join(records)}"
        )

    if len(records) < 2:
        raise SystemExit("FASTA must contain at least two records: original and design.")
    record_names = list(records)
    original_name = args.original_name or record_names[0]
    degenerate_name = args.degenerate_name or record_names[1]
    if original_name not in records:
        raise SystemExit(
            f"Original record {original_name!r} not found. Available: {', '.join(records)}"
        )
    if degenerate_name not in records:
        raise SystemExit(
            f"Design record {degenerate_name!r} not found. Available: {', '.join(records)}"
        )

    ori_full = records[original_name]
    deg_full = records[degenerate_name]
    if len(ori_full) != len(deg_full):
        raise SystemExit(
            f"Reference lengths differ: original={len(ori_full)}, design={len(deg_full)}. "
            "Sequences must be positionally aligned and have equal lengths."
        )

    start = (args.region_start - 1) if args.region_start is not None else 0
    end = args.region_end if args.region_end is not None else len(deg_full)
    if start < 0 or end > len(deg_full) or start >= end:
        raise SystemExit("Analyzed region must be within the sequences and have positive length.")
    region = (start, end)
    ori_seq = ori_full[start:end]
    deg_seq = deg_full[start:end]
    ori_protein = translate(ori_seq) if args.translate else None

    variable_positions = [
        i for i, ch in enumerate(deg_seq) if ch not in "ACGT"
    ]
    if not variable_positions:
        raise SystemExit("No IUPAC-degenerate positions found in the analyzed region.")

    mapping_seq = resolve_for_mapping(
        deg_full, ori_full, region
    )
    map_ref = outdir / "mapping_reference.fasta"
    with open(map_ref, "w") as fh:
        fh.write(">library_mapping_reference\n")
        for i in range(0, len(mapping_seq), 80):
            fh.write(mapping_seq[i:i+80] + "\n")

    bam = outdir / "aligned.sorted.bam"
    if not args.skip_mapping:
        run_mapping(map_ref, args.fastq, bam, args.preset)
    elif not bam.exists():
        raise SystemExit(f"--skip-mapping set, but {bam} does not exist.")

    print("[3/3] Calling library haplotypes...", file=sys.stderr)

    # Per-read categories and counters.
    total_primary = 0
    low_mapq = 0
    incomplete = 0
    invalid_design = 0
    assigned = 0

    hap_counts = Counter()
    invalid_counts = Counter()
    assignments = []

    # Per-designed-position counts among reads with a callable base at that site.
    pos_counts = {p: Counter() for p in variable_positions}
    pos_lowq = Counter()
    pos_deleted = Counter()
    pos_missing = Counter()

    # Fixed-position mismatch QC inside the analyzed region.
    fixed_positions = [i for i, ch in enumerate(deg_seq) if ch in "ACGT"]
    fixed_total = Counter()
    fixed_mismatch = Counter()
    fixed_alt = {i: Counter() for i in fixed_positions}

    seen_primary = set()

    for aln in sam_alignments_from_bam(bam):
        flag = aln["flag"]

        # unmapped, secondary, supplementary
        if flag & 0x4 or flag & 0x100 or flag & 0x800:
            continue
        if aln["qname"] in seen_primary:
            continue
        seen_primary.add(aln["qname"])
        total_primary += 1

        if aln["mapq"] < args.min_mapq:
            low_mapq += 1
            assignments.append((aln["qname"], "low_mapq", "", ""))
            continue

        obs = cigar_pairs(
            aln["pos"], aln["cigar"], aln["seq"], aln["qual"]
        )

        # Fixed-position mismatch QC.
        for p in fixed_positions:
            refpos = start + p
            if refpos not in obs:
                continue
            base, q = obs[refpos]
            if base == "-" or q < args.min_fixed_baseq:
                continue
            fixed_total[p] += 1
            if base != deg_seq[p]:
                fixed_mismatch[p] += 1
                fixed_alt[p][base] += 1

        called = []
        valid = True
        complete = True
        details = []

        for p in variable_positions:
            refpos = start + p
            code = deg_seq[p]
            allowed = IUPAC[code]

            if refpos not in obs:
                pos_missing[p] += 1
                complete = False
                details.append(f"{p+1}:missing")
                break

            base, q = obs[refpos]

            if base == "-":
                pos_deleted[p] += 1
                complete = False
                details.append(f"{p+1}:deletion")
                break

            if q < args.min_baseq:
                pos_lowq[p] += 1
                complete = False
                details.append(f"{p+1}:Q{q}")
                break

            pos_counts[p][base] += 1
            called.append(base)

            if base not in allowed:
                valid = False

        if not complete:
            incomplete += 1
            assignments.append(
                (aln["qname"], "incomplete_or_low_baseq", "", ";".join(details))
            )
            continue

        hap = "".join(called)

        if not valid:
            invalid_design += 1
            invalid_counts[hap] += 1
            assignments.append((aln["qname"], "invalid_design", hap, ""))
            continue

        assigned += 1
        hap_counts[hap] += 1
        assignments.append((aln["qname"], "assigned", hap, ""))

    # Build reconstructed DNA/protein for each valid haplotype.
    variant_rows = []
    protein_counts = Counter()
    protein_to_nt = defaultdict(list)

    def reconstruct(hap):
        chars = list(deg_seq)
        for p, b in zip(variable_positions, hap):
            chars[p] = b
        seq = "".join(chars)
        if any(ch not in "ACGT" for ch in seq):
            raise ValueError("Reconstructed sequence still contains ambiguity.")
        return seq

    for hap, n in hap_counts.most_common():
        dna = reconstruct(hap)
        protein = translate(dna) if args.translate else ""
        n_changes = nt_changes(ori_seq, dna)
        a_changes = aa_changes(ori_protein, protein) if args.translate else []
        if args.translate:
            protein_counts[protein] += n
            protein_to_nt[protein].append(hap)

        variant_rows.append({
            "reads": n,
            "frequency": (n / assigned) if assigned else 0.0,
            "haplotype": hap,
            "dna_sequence": dna,
            "nucleotide_change_count_vs_original": len(n_changes),
            "nucleotide_changes_vs_original": ";".join(n_changes),
            "protein_sequence": protein if args.translate else "",
            "aa_change_count_vs_original": len(a_changes),
            "amino_acid_changes_vs_original": ";".join(a_changes),
        })

    # variants.csv
    with open(outdir / "variants.csv", "w", newline="") as fh:
        fields = [
            "rank", "reads", "frequency", "haplotype", "dna_sequence",
            "nucleotide_change_count_vs_original",
            "nucleotide_changes_vs_original",
            "protein_sequence",
            "aa_change_count_vs_original",
            "amino_acid_changes_vs_original",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for rank, row in enumerate(variant_rows, start=1):
            w.writerow({"rank": rank, **row})

    # protein_variants.csv
    protein_rows = []
    for protein, n in protein_counts.most_common():
        aa_ch = aa_changes(ori_protein, protein)
        protein_rows.append({
            "reads": n,
            "frequency": (n / assigned) if assigned else 0.0,
            "protein_sequence": protein,
            "aa_change_count_vs_original": len(aa_ch),
            "amino_acid_changes_vs_original": ";".join(aa_ch),
            "number_of_distinct_nt_haplotypes": len(protein_to_nt[protein]),
            "nt_haplotypes": ";".join(protein_to_nt[protein]),
        })

    with open(outdir / "protein_variants.csv", "w", newline="") as fh:
        fields = [
            "rank", "reads", "frequency", "protein_sequence",
            "aa_change_count_vs_original", "amino_acid_changes_vs_original",
            "number_of_distinct_nt_haplotypes", "nt_haplotypes",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for rank, row in enumerate(protein_rows, start=1):
            w.writerow({"rank": rank, **row})

    # position_frequencies.csv
    with open(outdir / "position_frequencies.csv", "w", newline="") as fh:
        fields = [
            "amplicon_position_1based", "region_position_1based",
            "original_base", "design_code", "allowed_bases",
            "A_count", "C_count", "G_count", "T_count",
            "A_freq", "C_freq", "G_freq", "T_freq",
            "callable_reads", "low_quality_calls",
            "deletion_calls", "missing_calls",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for p in variable_positions:
            counts = pos_counts[p]
            callable_n = sum(counts[b] for b in "ACGT")
            row = {
                "amplicon_position_1based": start + p + 1,
                "region_position_1based": p + 1,
                "original_base": ori_seq[p],
                "design_code": deg_seq[p],
                "allowed_bases": "/".join(sorted(IUPAC[deg_seq[p]])),
                "callable_reads": callable_n,
                "low_quality_calls": pos_lowq[p],
                "deletion_calls": pos_deleted[p],
                "missing_calls": pos_missing[p],
            }
            for b in "ACGT":
                row[f"{b}_count"] = counts[b]
                row[f"{b}_freq"] = (counts[b] / callable_n) if callable_n else 0.0
            w.writerow(row)

    # fixed_position_mismatches.csv
    with open(outdir / "fixed_position_mismatches.csv", "w", newline="") as fh:
        fields = [
            "amplicon_position_1based", "region_position_1based",
            "reference_base", "callable_reads", "mismatch_reads",
            "mismatch_fraction", "A_alt", "C_alt", "G_alt", "T_alt",
        ]
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for p in fixed_positions:
            total = fixed_total[p]
            mism = fixed_mismatch[p]
            if total == 0:
                continue
            w.writerow({
                "amplicon_position_1based": start + p + 1,
                "region_position_1based": p + 1,
                "reference_base": deg_seq[p],
                "callable_reads": total,
                "mismatch_reads": mism,
                "mismatch_fraction": mism / total,
                "A_alt": fixed_alt[p]["A"],
                "C_alt": fixed_alt[p]["C"],
                "G_alt": fixed_alt[p]["G"],
                "T_alt": fixed_alt[p]["T"],
            })

    # invalid_haplotypes.csv
    with open(outdir / "invalid_haplotypes.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["haplotype", "reads"])
        for hap, n in invalid_counts.most_common():
            w.writerow([hap, n])

    # read_assignments.tsv
    with open(outdir / "read_assignments.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["read_name", "status", "haplotype", "details"])
        w.writerows(assignments)

    theoretical_size = 1
    for p in variable_positions:
        theoretical_size *= len(IUPAC[deg_seq[p]])

    # qc_summary.txt
    summary_lines = [
        "DEGENERATE LIBRARY VARIANT ANALYSIS",
        "===================================",
        f"Original FASTA record: {original_name}",
        f"Design FASTA record: {degenerate_name}",
        f"Original full reference length: {len(ori_full)} bp",
        f"Degenerate full reference length: {len(deg_full)} bp",
        f"Analyzed region: {start + 1}-{end} ({len(deg_seq)} bp)",
        f"Translation enabled: {args.translate}",
        f"Designed degenerate positions: {len(variable_positions)}",
        f"Theoretical nucleotide haplotypes from IUPAC design: {theoretical_size}",
        "",
        f"Primary mapped reads examined: {total_primary}",
        f"Reads below MAPQ {args.min_mapq}: {low_mapq}",
        f"Reads incomplete/low-quality at >=1 designed position: {incomplete}",
        f"Reads violating IUPAC design: {invalid_design}",
        f"Confidently assigned reads: {assigned}",
        "",
        f"Unique valid nucleotide haplotypes observed: {len(hap_counts)}",
        f"Unique protein variants observed: {len(protein_counts) if args.translate else 'not calculated'}",
        (
            "Fraction of confident reads represented by most abundant nucleotide variant: "
            f"{(hap_counts.most_common(1)[0][1] / assigned):.6f}"
            if assigned and hap_counts else
            "Fraction of confident reads represented by most abundant nucleotide variant: NA"
        ),
        "",
        "Interpretation note:",
        "Variant identity is determined ONLY by the designed IUPAC-degenerate positions.",
        "This prevents random sequencing errors at invariant positions from creating false variants.",
        "Unexpected mismatches at fixed analyzed positions are reported separately for QC.",
    ]
    with open(outdir / "qc_summary.txt", "w") as fh:
        fh.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines), file=sys.stderr)
    print(f"\nOutputs written to: {outdir.resolve()}", file=sys.stderr)


if __name__ == "__main__":
    main()
