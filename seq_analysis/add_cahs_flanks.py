#!/usr/bin/env python3
"""Add the supplied CAHS constant regions to a motif-only CSV (no dependencies).

Usage: python add_cahs_flanks.py input.csv output.csv
Preserves headers, row order, counts and any other columns.
"""
import argparse
import csv
from pathlib import Path

PREFIX = 'TEA'
SUFFIX = 'HLRDVEFRKDIVEMAIENQKKMIDVESRYAKKDMDRERVKV'
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--column', default='aa_variant', help='Sequence column name')
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        parser.error('Use a different output path to preserve the input.')

    # Read as strings so counts and other fields are not reformatted.
    with args.input.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames
        if not headers or args.column not in headers:
            parser.error(f'Missing sequence column: {args.column}')
        rows = list(reader)

    # Validate every row before writing; reject already-expanded sequences.
    for line, row in enumerate(rows, start=2):
        if None in row or any(value is None for value in row.values()):
            parser.error(f'Row {line}: malformed CSV row.')
        motif = row[args.column].strip().upper()
        if not motif or set(motif) - AMINO_ACIDS:
            parser.error(f'Row {line}: empty sequence or non-standard amino-acid characters.')
        if motif.startswith(PREFIX) and motif.endswith(SUFFIX):
            parser.error(f'Row {line}: constant regions appear to be present already.')
        # Concatenation adds no spaces, brackets or separators.
        row[args.column] = PREFIX + motif + SUFFIX

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    print(f'Wrote {len(rows)} sequences to {args.output}')


if __name__ == '__main__':
    main()
