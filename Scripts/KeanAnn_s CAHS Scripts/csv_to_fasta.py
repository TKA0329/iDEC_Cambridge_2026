"""
CSV to FASTA converter
Usage: python csv_to_fasta.py sequences.csv
       python csv_to_fasta.py sequences.csv --out output.fasta --prefix SEQ
"""

import sys
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert a sequence CSV to FASTA")
    parser.add_argument("input",           help="Input CSV file")
    parser.add_argument("--out",           default=None, help="Output FASTA path")
    parser.add_argument("--col",           default=None, help="Sequence column name (auto-detected if omitted)")
    parser.add_argument("--prefix",        default="SEQ", help="ID prefix for each entry (default: SEQ)")
    parser.add_argument("--id-col",        default=None, help="Optional column to use as FASTA IDs instead of numbering")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input)
    except Exception as e:
        sys.exit(f"Could not read file: {e}")

    # Auto-detect sequence column
    if args.col:
        seq_col = args.col
    else:
        candidates = [c for c in df.columns
                      if c.strip().lower() in
                      ("sequence", "seq", "protein", "aa_sequence", "protein_sequence")]
        if not candidates:
            sys.exit(f"Could not detect sequence column. Available: {list(df.columns)}\nUse --col to specify.")
        seq_col = candidates[0]
        print(f"Using column: '{seq_col}'")

    out_path = args.out or args.input.replace(".csv", ".fasta")

    with open(out_path, "w") as f:
        for i, row in enumerate(df[seq_col].astype(str), 1):
            seq = row.strip().upper().replace(" ", "")
            if not seq:
                continue
            # Use custom ID column if provided, otherwise number them
            if args.id_col and args.id_col in df.columns:
                seq_id = str(df[args.id_col].iloc[i - 1]).strip()
            else:
                seq_id = f"{args.prefix}_{i:04d}"
            f.write(f">{seq_id}\n{seq}\n")

    print(f"Done. {i} sequences written to: {out_path}")


if __name__ == "__main__":
    main()