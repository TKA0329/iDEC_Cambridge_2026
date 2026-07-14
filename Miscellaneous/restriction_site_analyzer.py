import itertools
import random
from typing import List, Dict, Generator

# IUPAC Degenerate base mapping
IUPAC: Dict[str, List[str]] = {
    'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T']
}

def rev_comp(seq: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return "".join(comp.get(base, base) for base in reversed(seq))

def total_combinations(seq: str) -> int:
    total = 1
    for base in seq:
        total *= len(IUPAC.get(base, [base]))
    return total

def generate_sequences(bases_options: List[List[str]]) -> Generator[str, None, None]:
    for p in itertools.product(*bases_options):
        yield "".join(p)

def random_seq(bases_options: List[List[str]]) -> str:
    return "".join(random.choice(opts) for opts in bases_options)

def format_number(n: int) -> str:
    if n >= 1e9: return f"{n/1e9:.2f}B"
    if n >= 1e6: return f"{n/1e6:.2f}M"
    if n >= 1e3: return f"{n/1e3:.1f}K"
    return f"{n:,}"

def main():
    print("=== Restriction Site Library Analyser ===")
    
    # 1. Get Degenerate Sequence Input
    degen_seq = input("Enter Degenerate codon sequence (5'→3'): ").strip().upper()
    if not degen_seq:
        print("Error: Sequence cannot be empty.")
        return
    invalid_degen = set(degen_seq) - set(IUPAC.keys())
    if invalid_degen:
        print(f"Error: Invalid IUPAC characters found: {', '.join(invalid_degen)}")
        return

    # 2. Get Restriction Enzyme Site Input
    enzyme_site = input("Enter Restriction enzyme site (5'→3'): ").strip().upper().replace(" ", "")
    if not enzyme_site:
        print("Error: Restriction site cannot be empty.")
        return
    invalid_enzyme = set(enzyme_site) - {'A', 'C', 'G', 'T'}
    if invalid_enzyme:
        print(f"Error: Restriction site must use standard bases (A, C, G, T). Got: {', '.join(invalid_enzyme)}")
        return
        
    if len(enzyme_site) > len(degen_seq):
        print("Error: Restriction site is longer than your degenerate sequence.")
        return

    # 3. Get Reverse Complement Option
    rc_input = input("Also check reverse complement? (Y/n): ").strip().lower()
    check_rc = rc_input not in ['n', 'no']

    # --- Analysis setup ---
    bases_options = [IUPAC[b] for b in degen_seq]
    site_rc = rev_comp(enzyme_site)
    total_combos = total_combinations(degen_seq)
    
    MAX_ENUM = 500000
    SAMPLE_N = 150000
    
    def contains_site(s: str) -> bool:
        if enzyme_site in s: return True
        if check_rc and enzyme_site != site_rc and site_rc in s: return True
        return False

    match_count = 0
    is_sampled = False

    print("\nAnalysing library... please wait...")
    
    if total_combos <= MAX_ENUM:
        for s in generate_sequences(bases_options):
            if contains_site(s):
                match_count += 1
    else:
        is_sampled = True
        for _ in range(SAMPLE_N):
            if contains_site(random_seq(bases_options)):
                match_count += 1
        match_count = round((match_count / SAMPLE_N) * total_combos)

    proportion = match_count / total_combos
    percentage = proportion * 100

    # --- Output Results ---
    print("\n" + "="*40)
    print("           ANALYSIS RESULTS           ")
    print("="*40)
    print(f"Method used:        {'~ESTIMATED (Sampled)' if is_sampled else 'EXACT (Enumerated)'}")
    print(f"Total Library Size: {format_number(total_combos)}")
    print(f"Sequences w/ site:  {format_number(match_count)}")
    print(f"Frequency:          {percentage:.2f}%")
    print("-"*40)
    print(f"Enzyme site target: 5'-{enzyme_site}-3'")
    if check_rc and enzyme_site != site_rc:
        print(f"RC site checked:    5'-{site_rc}-3'")
    elif enzyme_site == site_rc:
        print("Palindromic site — Reverse complement is identical")
    print("="*40 + "\n")

if __name__ == "__main__":
    main()