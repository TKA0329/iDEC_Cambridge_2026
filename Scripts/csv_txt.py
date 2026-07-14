import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Convert a sequence CSV to FASTA format.')
parser.add_argument('input_file', nargs='?', default='sequences.csv',
                     help='Path to input CSV file (default: sequences.csv)')
parser.add_argument('output_file', nargs='?', default='output.txt',
                     help='Path to output FASTA file (default: output.txt)')
args = parser.parse_args()

df = pd.read_csv(args.input_file)
with open(args.output_file, 'w') as f:
    for i, seq in enumerate(df['sequence'], 1):
        f.write(f'>{i}\n{seq}\n')