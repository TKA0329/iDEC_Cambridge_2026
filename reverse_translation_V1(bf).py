from itertools import product
import sys

# =========================================================
# IUPAC ambiguity codes
# =========================================================

IUPAC = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

# =========================================================
# Standard genetic code
# =========================================================

CODON_TABLE = {

    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',

    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',

    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',

    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# =========================================================
# FASTA parser
# =========================================================

def read_fasta(filename):

    sequences = []
    current = []

    with open(filename) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):

                if current:
                    sequences.append(''.join(current))
                    current = []

            else:
                current.append(line)

        if current:
            sequences.append(''.join(current))

    return sequences

# =========================================================
# Expand degenerate codon
# =========================================================

def expand_degenerate_codon(codon):

    possibilities = []

    for base in codon:
        possibilities.append(IUPAC[base])

    concrete_codons = []

    for combo in product(*possibilities):
        concrete_codons.append(''.join(combo))

    return concrete_codons

# =========================================================
# Determine amino acids encoded
# =========================================================

def codon_to_amino_acids(degenerate_codon):

    concrete_codons = expand_degenerate_codon(degenerate_codon)

    amino_acids = set()

    for codon in concrete_codons:

        aa = CODON_TABLE[codon]

        amino_acids.add(aa)

    return amino_acids

# =========================================================
# Generate all possible degenerate codons
# =========================================================

def generate_all_degenerate_codons():

    symbols = list(IUPAC.keys())

    codons = []

    for a in symbols:
        for b in symbols:
            for c in symbols:

                codons.append(a + b + c)

    return codons

# =========================================================
# Build lookup table
# =========================================================

def build_lookup_table():

    lookup = {}

    all_codons = generate_all_degenerate_codons()

    for codon in all_codons:

        aas = codon_to_amino_acids(codon)

        lookup[codon] = aas

    return lookup

# =========================================================
# Scoring function
# =========================================================

def score_candidate(target_aas, encoded_aas):

    missing = len(target_aas - encoded_aas)

    unwanted = len(encoded_aas - target_aas - {'*'})

    stops = 1 if '*' in encoded_aas else 0

    size_penalty = len(encoded_aas)

    score = (
        missing * 1000 +
        unwanted * 10 +
        stops * 100 +
        size_penalty
    )

    return score

# =========================================================
# Find best codon
# =========================================================

def find_best_codon(target_aas, lookup_table):

    best_codon = None
    best_score = float('inf')
    best_encoded = None

    for codon, encoded_aas in lookup_table.items():

        if not target_aas.issubset(encoded_aas):
            continue

        score = score_candidate(target_aas, encoded_aas)

        if score < best_score:

            best_score = score
            best_codon = codon
            best_encoded = encoded_aas

    return best_codon, best_encoded, best_score

# =========================================================
# Analyze alignment
# =========================================================

def analyze_alignment(sequences):

    if len(sequences) == 0:
        raise ValueError("No sequences found in FASTA file.")

    length = len(sequences[0])

    for seq in sequences:

        if len(seq) != length:
            raise ValueError(
                "Sequences are not aligned.\n"
                "All sequences must have the same length."
            )

    columns = []

    for i in range(length):

        aas = set()

        for seq in sequences:

            aa = seq[i]

            if aa != '-':
                aas.add(aa)

        columns.append(aas)

    return columns

# =========================================================
# Main
# =========================================================

def main():

    if len(sys.argv) < 2:

        print("\nUsage:")
        print("python reverse_translation.py input.fasta\n")

        return

    fasta_file = sys.argv[1]

    print(f"\nReading FASTA file: {fasta_file}")

    sequences = read_fasta(fasta_file)

    print(f"Loaded {len(sequences)} sequences")

    lookup = build_lookup_table()

    print(f"Generated {len(lookup)} degenerate codons")

    columns = analyze_alignment(sequences)

    print("\nOptimized Degenerate Codons:\n")

    final_sequence = ""

    intended_diversity = 1
    actual_diversity = 1

    for i, target_aas in enumerate(columns):

        best_codon, encoded, score = find_best_codon(
            target_aas,
            lookup
        )

        final_sequence += best_codon

        intended_diversity *= len(target_aas)

        actual_diversity *= len(encoded - {'*'})

        print(f"Position {i+1}")
        print(f"  Target amino acids : {sorted(target_aas)}")
        print(f"  Best codon         : {best_codon}")
        print(f"  Encoded amino acids: {sorted(encoded)}")
        print()

    off_target = actual_diversity - intended_diversity

    efficiency = (
        intended_diversity / actual_diversity
    ) * 100

    nnk_diversity = 20 ** len(columns)

    nnk_off_target = nnk_diversity - intended_diversity

    nnk_efficiency = (
        intended_diversity / nnk_diversity
    ) * 100

    print("=================================================")
    print("FINAL DEGENERATE DNA SEQUENCE")
    print("=================================================\n")

    print(final_sequence)
    print()

    print("=================================================")
    print("OPTIMIZED LIBRARY STATISTICS")
    print("=================================================\n")

    print(f"Intended protein variants : {intended_diversity}")
    print(f"Encoded protein variants  : {actual_diversity}")
    print(f"Off-target variants       : {off_target}")
    print(f"Library efficiency        : {efficiency:.6f}%")

    print("\n=================================================")
    print("NNK COMPARISON")
    print("=================================================\n")

    print(f"NNK protein diversity     : {nnk_diversity}")
    print(f"NNK off-target variants   : {nnk_off_target}")
    print(f"NNK efficiency            : {nnk_efficiency:.12f}%")

# =========================================================

if __name__ == '__main__':
    main()