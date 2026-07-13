import csv
import sys

def fasta_to_csv(fasta_file, csv_file):
    sequences = []
    with open(fasta_file) as f:
        seq = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                seq = ""
            else:
                seq += line
        if seq:
            sequences.append(seq)

    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence"])
        for seq in sequences:
            writer.writerow([seq])

    print(f"Written {len(sequences)} sequences to {csv_file}")

if __name__ == "__main__":
    fasta_to_csv(sys.argv[1], sys.argv[2])