"""
workflow - code's logic 

1. parse fastq, get sequences

2. locate frame-defining anchor -> using edlib instead of .find() because 
sequencing is fuzzy and want to tolerate some errors in string matching 

3. since only the 19AA motif 1 is variable (with degeneracy), will use the rest of the constant region
to determine if its a "good" sequence; if its full of indels/weird stuff, discard the whole sequence

4. pull a fixed length 57 nt (19-codon) window from that frame definining anchor from the "good"
sequences. 

Evaluate the base qualities in the 57 nt and flag 

5. Translate these survivors to amino acids, count with a counter using @ dataclass

6. Error-correct: greedily merge a low count variant into a more abundant neighbour.
Determined by a hamming distance (used to determine no. of mismatches) < 2
Why? because an error is far more likely than a second independent clone landing 1-2 residues 
from the dominant one in a motif this degenerate.

8. Report the top cluusters, will prioritise sequences with count >= 8 

# implementing two checks:
1. cross check with degenerate sequence and check if the regions that are not supposed to be 
mutated are mutated. If they are mutated -> discard; if not -> keep 
2. quality score weighted analysis -> ?? 

usage 
python sequencing_data_analyzerV2.py VRZQTP_1_sample_1.fastq cahs_dna_variants.csv cahs_aa_variants.csv valid_var.txt discarded_var_at_find_anchor.csv

"""

import edlib 
from dataclasses import dataclass, field
from tqdm import tqdm
import argparse 
from Bio.Seq import Seq 
from collections import Counter 
import csv

# HYPERPARAMS
LEFT_ANCHOR = "TAACAGGAGGAATTAACCATGACAGAAGCA" # 30 nt 
CAHS_SEQ_WHOLE = "GTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGACAGAAGCAAKRMKGARAMRAATGRVAARAGMAADGRRSARARTAAGAAAASAAWTGRMAAAAATGCACCTAAGAGACGTAGAATTCAGAAAAGACATAGTAGAAATGGCAATAGAAAACCAAAAAAAAATGATAGACGTAGAAAGCAGATACGCAAAAAAAGACATGGACAGAGAAAGAGTAAAAGTATAAAGCTCGAGCGAAGCTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGTCTCCAGCTTGGCTGTTTTGGCGGATGA"
DEGEN_SEQ = "AKRMKGARAMRAATGRVAARAGMAADGRRSARARTAAGAAAASAAWTGRMAAAAATG"
RIGHT_ANCHOR_CHECK = "CACCTAAGAGACGTAGAATTCAGAAAAGAC" # 30 nt 

MOTIF_LEN_AA = 19
MOTIF_LEN_NT = MOTIF_LEN_AA * 3 

FUZZY_K = 8 # allows 8 differences between LEFT_ANCHOR and what appears in read 
QC_TOLERANCE_K = 8 # allows 8 differences between LEFT_ANCHOR and what appears in read
CLUSTER_RATIO = 3 # a variant is absorbed if the neighbour >= this many times its count 
# rationale: if a neighbour has 20 reads, and a similar sequence has 1 read, 20/1 
# -> very likely to be sequencing error
# but not fool proof, so has to satisfy this condition:
CLUSTER_MAX_DIST = 2 # within this many AA substitutions of it 

IUPAC_DEGEN_BASES = {
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


def parse_fastq(filepath):
    """
    A generator that yields (read_id, sequence and quality_scores) tuple
    """
    with open(filepath, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            sequence = f.readline().rstrip("\n")
            plus = f.readline()
            quality = f.readline().rstrip("\n")
            header = header.rstrip("\n")
        

            if not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ: expected '@', got: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ: expected '+' separator, got {plus!r}")
            if len(sequence) != len(quality):
                raise ValueError(f"Malformed FASTQ: sequence/quality length mismatch for read {header!r}")
            read_id = header[1:] # strip the leading "@"

            yield read_id, sequence, quality #using yield to keep giving values as a loop progresses

# reverse_complement
def reverse_complement(sequence):
    complement = {"A":"T",
                  "T":"A",
                  "C":"G",
                  "G":"C",
                  "N": "N"}
    return "".join(complement[base] for base in reversed(sequence))

@dataclass
class Metrics:
    total_reads: int = 0 # total number of reads provided by the sequencing company 
    no_left_anchor: int = 0 
    short_window: int = 0
    indel_in_right_window: int = 0
    stop_codon: int = 0 
    called: int = 0 
    aa_counts: Counter = field(default_factory=Counter)
    mismatch_with_degenerate_seq: int=0


metrics = Metrics()
def find_anchor(fastq_path):
    """
    make sure that the sequence passes a series of checks
    then only motif 1 of that sequence gets translated 
    """
    
    rows = []
    list_nt_variant = []
    list_invalid=[]
    for read_id, seq, qual in tqdm(parse_fastq(fastq_path)):
        metrics.total_reads += 1
        strand = "+"
        
        # locating where LEFT ANCHOR is, and from the first nt after LEFT ANCHOR ended,
        # identify motif 1, reverse_complementing strands wehre necessary
        la = edlib.align(LEFT_ANCHOR, seq, mode="HW",
                         task="locations", k=FUZZY_K) # allowing 8 differences, semi-global alignment with HW

        if la["editDistance"] == -1:
            seq = reverse_complement(seq)
            strand = "-"
            la = edlib.align(LEFT_ANCHOR, seq, mode="HW",
                             task="locations", k = FUZZY_K)
            if la["editDistance"] == -1:
                metrics.no_left_anchor += 1 
                list_invalid.append({"read_id": read_id,
                                     "seq": seq,
                                     "reason": "no_left_anchor"})
                continue 
        frame_start = la["locations"][0][1] + 1 # because edlib's locations are (start, end)
        # inclusive of the last matched base. 
        nt_variant = seq[frame_start:frame_start+MOTIF_LEN_NT]
        
        if len(nt_variant) < MOTIF_LEN_NT:
            metrics.short_window += 1
            list_invalid.append({"read_id":read_id,
                                 "seq":seq,
                                 "reason":"short_window"}) 
            
            
            continue

        # does the constant post_motif sequence sit exactly where a clean, 57 nt, 
        # no indels window predicts ? (+/- a few nt diff)

        downstream_window = seq[frame_start+MOTIF_LEN_NT:frame_start+MOTIF_LEN_NT+len(RIGHT_ANCHOR_CHECK)] 
        # not sure if it should have more freedom of letting it slide up and down?
        # but even if it does land 1 or 2 nucleotides further upstream of the variant
        # we wouldn't know where the insertion or deletion takes place and its impossible to guess
        # so sequence should still be discarded. 
        ra = edlib.align(RIGHT_ANCHOR_CHECK, downstream_window, mode="HW", task="locations",
                         k=QC_TOLERANCE_K)
        
        if ra["editDistance"] == -1:
            metrics.indel_in_right_window += 1
            list_invalid.append({"read_id":read_id,
                                 "seq":seq,
                                 "reason": "indel_in_right_window"})
            continue

        list_nt_variant.append(nt_variant)

        rows.append({"strand": strand,
                            "nt_variant": nt_variant})

        # shouldnt use zip for list_invalid because it'll cap how many rows you get based on which
        # is rarest and will lose the vast majority of no_left_anchor / short windows etc. 
    
    return metrics, rows, list_nt_variant, list_invalid

    

def discard_invalid_variant(list_nt_variant: list, degen_seq=DEGEN_SEQ):
    """
    discard sequences with mutation regions that are not supposed to be mutated
    by comparing them with the degenerate sequence and the bases allowed 
    at each degenerate position
    """
    kept_list = []
    discarded_list = []
    for variant in list_nt_variant:
        is_valid = True
        if len(variant) != len(degen_seq):
            discarded_list.append(variant)
            continue
        for obs_base, degen_code in zip(variant, degen_seq):
            allowed = IUPAC_DEGEN_BASES.get(degen_code.upper())
            if allowed is None or obs_base.upper() not in allowed:
                is_valid = False
                break
        if is_valid:
            kept_list.append(variant)
        else:
            discarded_list.append(variant)
            metrics.mismatch_with_degenerate_seq += 1
    return kept_list, discarded_list, metrics

def convert_to_aa(list_nt_variant):
    aa_variants = []
    for nt_variant in list_nt_variant:
        
        aa = str(Seq(nt_variant).translate())
        if "*" in aa:
            metrics.stop_codon += 1
            continue 

        metrics.called += 1 
        metrics.aa_counts[aa] += 1
        aa_variants.append({"aa": aa})
    
    return aa_variants


def hamming(a: str, b:str) -> int:
    return sum(x != y for x,y in zip(a,b))

def error_correct_clusters(aa_counts: Counter, 
                           ratio = CLUSTER_RATIO, max_dist=CLUSTER_MAX_DIST):
    """
    Greedy directional clustering
    """
    items = sorted(aa_counts.items(), key=lambda x: -x[1]) 
    # sort the dict from largest to smallest
    # also turn dict into a list-like structure e.g. {"A":"B"} -> [("A":"B")]
    
    
    clusters = [] # shud look like: (representative, members=[(seq, count), ...])
    for seq, count in items: 
        for i, (rep, members) in enumerate(clusters):
            total_count = sum(c for _, c in members) # calculates total no. of reads of that cluster
            if total_count >= ratio * count and hamming(seq, rep) <= max_dist: 
                # gets appended to the same group if it doesn't pass the check
                # i.e. its total_count/count >= ratio and hamming distance <= 2 
                members.append((seq, count)) 
                break
        else: # otherwise if its diff enough, gets to be a new rep 
            clusters.append((seq, [(seq, count)]))
        
       

    clusters.sort(key = lambda c:-sum(count for _,count in c[1]))
    return clusters 

def main():
    parser = argparse.ArgumentParser(
    description = "Convert a FASTQ file into a csv.")
    parser.add_argument("input_fastq", help="Path to the input .fastq file")
    parser.add_argument("dna_var_csv", 
                        help="Path to write the output .csv file containing the dna variants that pass the find_anchor check")
    parser.add_argument("aa_var_csv", 
                        help="Path to write the output .csv file containing the aa acid variants before clustering ")
    parser.add_argument("valid_var_txt", 
                        help="Path to write the output txt file containing the valid sequences after all the checks + their frequency")
    parser.add_argument("invalid_var_discarded_at_find_anchor_csv", 
                        help="Path to write the output csv file containing the variants discarded at the find_anchor check")
    args = parser.parse_args()

    metrics, rows, list_variant, list_invalid = find_anchor(args.input_fastq)
    valid_variant, discarded_variant, metrics = discard_invalid_variant(list_variant)
    valid_aa = convert_to_aa(valid_variant)
    total = metrics.called 

    print("-----------------read-level QC----------------------")
    print(f"total reads: {metrics.total_reads}")
    print(f"no left anchor found: {metrics.no_left_anchor}")
    print(f"window too short: {metrics.short_window}")
    print(f"indel inside motif window: {metrics.indel_in_right_window}")
    print(f"premature stop codon: {metrics.stop_codon}")
    print(f"mismatch with the degenerate seq: {metrics.mismatch_with_degenerate_seq}")
    print(f"cleanly called: {metrics.called} ({metrics.called/metrics.total_reads* 100:.0f} %)")
    print(f"raw distinct AA strings: {len(metrics.aa_counts)}")
    clusters = error_correct_clusters(metrics.aa_counts)
    print(f"error-corrected clusters (Hamming <= {CLUSTER_MAX_DIST}, >= {CLUSTER_RATIO} x abundance ratio): {len(clusters)}")
    print(f"--------------- top 10 variants (error-corrected) -----------------")
    
    for rep, members in clusters[:10]: # must cover > 80 % of the sequences 
        count = sum(c for _, c in members) # counting no. of reads in the top 10 variants
        percentage = (count/total)*100 # percentage of x sequence that exists in the sequencing data
                                         
        print(f" {count:4d} reads ({percentage:5.1f}%) {rep}")
        print(f"-> {len(members)} raw variant(s) merged")
    with open(args.valid_var_txt, "w") as f:
        for rep, members in clusters: # all the sequences available
            count = sum(c for _, c in members) 
            
            percentage = (count/total)*100 # percentage of x sequence that exists in the sequencing data
            f.write(f"Sequence: {rep}; percentage out of valid variants {total}: {percentage:.1f} % \n")
    print(f"-------------------------------------------------------------")
    
    with open(args.dna_var_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["strand", "nt_variant"])
        w.writeheader()
        w.writerows(rows)

    
    with open(args.aa_var_csv, "w", newline="") as file:
        write = csv.DictWriter(file, fieldnames=["aa_variant", "count"])
        write.writeheader()
        for aa_seq, count in metrics.aa_counts.items():
            write.writerow({"aa_variant": aa_seq, "count":count})

    with open(args.invalid_var_discarded_at_find_anchor_csv, "w", newline="") as file_:
        write = csv.DictWriter(file_, fieldnames=["read_id", "seq", "reason"])
        write.writeheader()
        write.writerows(list_invalid) 

if __name__ == "__main__":
    main()
 




            


    



        
        









