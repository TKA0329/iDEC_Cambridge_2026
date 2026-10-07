# """
# usage 
# python sequencing_data_analyzerV3.py VRZQTP_4_sample_4.fastq ssb_dna_variants.csv ssb_aa_variants.csv valid_var.txt discarded_var_at_find_anchor.csv discarded_degen_mis.csv

"""
usage 
python sequencing_data_analyzerV3.py VRZQTP_1_sample_1.fastq cahs_dna_variants.csv cahs_aa_variants.csv valid_var.txt discarded_var_at_find_anchor.csv discarded_degen_mis.csv

"""

import edlib 
from dataclasses import dataclass, field, replace, asdict
from tqdm import tqdm
import argparse 
from Bio.Seq import Seq 
from collections import Counter 
import csv
import re






# HYPERPARAMS
LEFT_ANCHOR = "CTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGACAGAAGCA" # 60 nt 
CAHS_SEQ_WHOLE = "GTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTGGGCTAACAGGAGGAATTAACCATGACAGAAGCAAKRMKGARAMRAATGRVAARAGMAADGRRSARARTAAGAAAASAAWTGRMAAAAATGCACCTAAGAGACGTAGAATTCAGAAAAGACATAGTAGAAATGGCAATAGAAAACCAAAAAAAAATGATAGACGTAGAAAGCAGATACGCAAAAAAAGACATGGACAGAGAAAGAGTAAAAGTATAAAGCTCGAGCGAAGCTTGGGCCCGAACAAAAACTCATCTCAGAAGAGGATCTGAATAGCGCCGTCGACCATCATCATCATCATCATTGAGTTTAAACGGTCTCCAGCTTGGCTGTTTTGGCGGATGA"
DEGEN_SEQ = "AKRMKGARAMRAATGRVAARAGMAADGRRSARARTAAGAAAASAAWTGRMAAAAATG"
RIGHT_ANCHOR_CHECK = "CACCTAAGAGACGTAGAATTCAGAAAAGACATAGTAGAAATGGCAATAGAAAACCAAAAA" # 60 nt 

MOTIF_LEN_AA = 19
MOTIF_LEN_NT = MOTIF_LEN_AA * 3 

FUZZY_K = 16 # allows 8 differences between LEFT_ANCHOR and what appears in read 
QC_TOLERANCE_K = 16 # allows 8 differences between LEFT_ANCHOR and what appears in read
CLUSTER_RATIO = 3 # a variant is absorbed if the neighbour >= this many times its count 
# rationale: if a neighbour has 20 reads, and a similar sequence has 1 read, 20/1 
# -> very likely to be sequencing error
# but not fool proof, so has to satisfy this condition:
CLUSTER_MAX_DIST = 2 # within this many AA substitutions of it 

RESCUE_K = FUZZY_K + QC_TOLERANCE_K

CONFIDENT_Q = 30
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
    confident_mismatch: int= 0
    ambiguous_low_q: int=0
    mismatch_corrected: int=0
    rescued_indels: int=0
    

metrics = Metrics()


# The point of Read is to give each sequencing read a standard container for all the 
# information you want to carry through your pipeline.
# Without Read, you'd probably have lots of separate variables/lists, and its 
# hard to make sure everything stays aligned. 
@dataclass
class Read:
    read_id: str
    strand: str
    nt: str # the motif window
    qual: str # quality string; same length as nt 
    source: str # "direct" or "rescued"

def build_iupac_equalities():
    pairs = []
    for code, bases in IUPAC_DEGEN_BASES.items():
        for base in bases:
            if base != code:
                pairs.append((code, base))
    return pairs

IUPAC_EQUALITIES = build_iupac_equalities()

REF_WINDOW = LEFT_ANCHOR + DEGEN_SEQ + RIGHT_ANCHOR_CHECK   # 30 + 57 + 30 = 117 nt
MOTIF_REF_START = len(LEFT_ANCHOR)                          # 30
MOTIF_REF_END   = MOTIF_REF_START + MOTIF_LEN_NT            # 87

def find_anchor(fastq_path):
    """
    make sure that the sequence passes a series of checks
    then only motif 1 of that sequence gets translated 
    """
    
    reads = []
    # list_nt_variant = []
    list_invalid=[]
    for read_id, seq, qual in tqdm(parse_fastq(fastq_path)):
        orig_seq, orig_qual = seq, qual
        metrics.total_reads += 1
        strand = "+"
        
        # locating where LEFT ANCHOR is, and from the first nt after LEFT ANCHOR ended,
        # identify motif 1, reverse_complementing strands wehre necessary
        la = edlib.align(LEFT_ANCHOR, seq, mode="HW",
                         task="locations", k=FUZZY_K) # allowing 8 differences, semi-global alignment with HW

        if la["editDistance"] == -1: # meaning it doesn't match 
            seq = reverse_complement(seq)
            qual = qual[::-1]
            strand = "-"
            la = edlib.align(LEFT_ANCHOR, seq, mode="HW",
                             task="locations", k = FUZZY_K)
            if la["editDistance"] == -1:
                metrics.no_left_anchor += 1 
                list_invalid.append({"read_id": read_id,
                                     "seq": orig_seq, # storing the rev comp strand 
                                     "reason": "no_left_anchor",
                                     "qual":orig_qual,
                                     "strand": "?"}) # problem with this: if its rev comp and doesnt match, will remain rev comp
                # affects downstream processes and so strand is seen as "?" instead of "+" or "-" and orig_seq/orig_qual
                continue 
        
        frame_start = la["locations"][0][1] + 1 # because edlib's locations are (start, end)
        # inclusive of the last matched base. 

        # if rev comp works then will be processed as rev comp strand from here onwards
        nt_variant = seq[frame_start:frame_start+MOTIF_LEN_NT]

        # if idp is too short 
        if len(nt_variant) < MOTIF_LEN_NT: # by this point we already know the correct orientation of the strand 
            metrics.short_window += 1
            list_invalid.append({"read_id":read_id,
                                 "seq":seq,
                                 "reason":"short_window",
                                 "qual": qual,
                                 "strand": strand}) 
            
            
            continue

        # does the constant post_motif sequence sit exactly where a clean, LEN_MOTIF_NT nt, 
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
                                 "reason": "indel_in_right_window",
                                 "qual": qual,
                                 "strand": strand})
            continue

        # list_nt_variant.append(nt_variant)

        # good reads get appended to rows 
        reads.append(Read(read_id, strand, nt_variant, qual[frame_start:frame_start+MOTIF_LEN_NT], "direct"))

        # shouldnt use zip for list_invalid because it'll cap how many rows you get based on which
        # is rarest and will lose the vast majority of no_left_anchor / short windows etc. 
    
    return metrics, reads, list_invalid

def rescue_indel_read(list_invalid, ref_window, motif_ref_start, motif_ref_end) -> list: 

    """
    REF_WINDOW:   [ LEFT_ANCHOR (30) ][   MOTIF (57)   ][ RIGHT_ANCHOR_CHECK (30) ]
    ref index:    0                   30               87                          117

    READ:   ...junk... [ left flank ][    motif    ][ right flank ] ...junk...
                       ^             ^              ^              ^
                locations[0][0]  motif_start_q  motif_end_q   locations[0][1]
    
    """
    rescued_list = []
    for discard in list_invalid: 
        # list_invalid is [{"read_id":read_id, "seq":seq, "reason": "indel_in_right_window", "qual": qual,"strand": strand}]

        # rescue_indel_read, replacing the single align call
        # "?" means unsure of orientation of strand 
        if discard["strand"] == "?":
            candidates = [("+", discard["seq"], discard["qual"]),
                        ("-", reverse_complement(discard["seq"]), discard["qual"][::-1])]
        else:
            candidates = [(discard["strand"], discard["seq"], discard["qual"])]

        # finds the best alignment for the sequence 
        best = None
        for strand, s, q in candidates:
            aln = edlib.align(ref_window, s, mode="HW", task="path",
                            k=RESCUE_K, additionalEqualities=IUPAC_EQUALITIES)
            # if mismatch or number of mismatches < current best 
            # if is -1, unrecoverable because mismatches more than RESCUE_K i.e. completely unalignable
            if aln["editDistance"] != -1 and (best is None or aln["editDistance"] < best[0]["editDistance"]):
                best = (aln, strand, s, q)

        if best is None:
            continue
        aln, strand, seq, qual = best
        
        # CIGAR walk as before, but slice seq and qual, and use strand in the Read
        
        ref_pos = 0 # index into ref_window
        read_pos = aln["locations"][0][0] # where actual correct read starts 
        motif_start_q = motif_end_q = None
        
        indel_at_motif = False 

        # regex 
        for length, op in re.findall(r"(\d+)([=XID])", aln["cigar"]):
            # e.g. suppose aln["cigar"] = "3=1I5=2D"
            # re.findall returns [("3", "="), ("1", "I"), ("5", "="),("2", "D")]

            for _ in range(int(length)):
                # it has reached the start of IDP 
                if ref_pos == motif_ref_start and motif_start_q is None:
                    # motif_ref_start is 30 nt
                    motif_start_q = read_pos # both read and ref pos move tgt as CIGAR walks 
                    # so aln["locations"][0][0] tells you overall where along the read motif_start is
                
                # it has reached the end of IDP 
                if ref_pos == motif_ref_end and motif_end_q is None:
                    # motif_ref_end is 30 + length of idp
                    motif_end_q = read_pos

                # if indels inside IDP 
                if op in "ID" and motif_ref_start <= ref_pos <= motif_ref_end:
                    indel_at_motif = True # boundary counts as ambiguous 


                if op in "=X":
                    ref_pos += 1; read_pos += 1

                elif op == "D": 
                    read_pos += 1

                elif op == "I":
                    ref_pos += 1 
        
        if indel_at_motif or motif_start_q is None or motif_end_q is None:
            continue

        motif = seq[motif_start_q:motif_end_q]

        if len(motif) == MOTIF_LEN_NT:
            
            rescued_list.append(Read(discard["read_id"], strand, 
                                     motif, qual[motif_start_q:motif_end_q], "rescued"))
            
            metrics.rescued_indels += 1
            
        else:
            continue
    return rescued_list

def phred(ch):
    return ord(ch) - 33

def discard_invalid_variant(reads, degen_seq=DEGEN_SEQ):
    """
    discard sequences with mutation regions that are not supposed to be mutated
    by comparing them with the degenerate sequence and the bases allowed 
    at each degenerate position
    """
    
    kept, discarded = [], []
    for r in reads:
        # checking the length of the reads that have passed the anchor checks 
        if len(r.nt) != len(degen_seq):
            discarded.append({"read_id": r.read_id, "dna_seq": r.nt, "reason": "length"})
            continue

        
        fixed, reason = [], None
        # check constant regions
        # obs - observed; code - correct base in consensus seq; q - quality 
        for obs, code, q in zip(r.nt, degen_seq, r.qual):
            # if the base is in allowed, append to the fixed list 
            allowed = IUPAC_DEGEN_BASES[code.upper()]
            if obs.upper() in allowed:
                fixed.append(obs)

            # if its not in allowed, if its phred socre is high (high confidence), break the loop 
            elif phred(q) >= CONFIDENT_Q:
                reason = "confident_mismatch"      # real deviation -> discard
                metrics.confident_mismatch += 1
                break
            # if its not in allowed, and its phred score < Q, and its a constant base, append the correct one
            elif len(allowed) == 1:
                fixed.append(next(iter(allowed)))  # low-Q error at a fixed position -> correct
            else:
                # because cannot deterministically know which base is correct when there are multiple you can choose from
                reason = "ambiguous_low_q"         # can't tell which allowed base -> discard
                metrics.ambiguous_low_q += 1
                break
        # if there's a reason, this means its discarded
        if reason:
            discarded.append({"read_id": r.read_id, "dna_seq": r.nt, "reason": reason})
            metrics.mismatch_with_degenerate_seq += 1
        else:
            kept.append(replace(r, nt="".join(fixed)))  
            metrics.mismatch_corrected += 1
    return kept, discarded

def convert_to_aa(list_nt_variant):
    aa_variants = []
    for r in list_nt_variant:
        
        aa = str(Seq(r.nt).translate())
        if "*" in aa:
            metrics.stop_codon += 1
            continue 

        metrics.called += 1 
        metrics.aa_counts[aa] += 1
        aa_variants.append({"read_id": r.read_id,
        "aa":aa})
    
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
    parser.add_argument("discarded_at_degen_mismatch", help="Path to write the output csv file containing the variants discarded at degen mismatch")
    args = parser.parse_args()

    metrics, reads, invalid = find_anchor(args.input_fastq)
    rescued = rescue_indel_read(invalid, REF_WINDOW, MOTIF_REF_START, MOTIF_REF_END)
    reads = reads + rescued

    kept, discarded = discard_invalid_variant(reads)
    valid_aa = convert_to_aa(kept)      # inside: Seq(r.nt).translate()
    
    total = metrics.called 
    

    print("-----------------read-level QC----------------------")
    print(f"total reads: {metrics.total_reads}")
    print(f"no left anchor found: {metrics.no_left_anchor}")
    print(f"window too short: {metrics.short_window}")
    print(f"indel inside motif window: {metrics.indel_in_right_window}")
    print(f"no. of seq rescued from rescue_indel_read: {metrics.rescued_indels}")
    print(f"premature stop codon: {metrics.stop_codon}")
    print(f"mismatch with the degenerate seq: {metrics.mismatch_with_degenerate_seq}")
    print(f"        mismatch due to high confidence in mismatch: {metrics.confident_mismatch}")
    print(f"        mismatch at degenerate regions (ambiguous low q): {metrics.ambiguous_low_q}")
    print(f"        mismatch corrected: {metrics.mismatch_corrected}")
    print(f"cleanly called: {metrics.called} ({metrics.called/metrics.total_reads* 100:.0f} %)")
    print(f"raw distinct AA strings (unclustered): {len(metrics.aa_counts)}")
    clusters = error_correct_clusters(metrics.aa_counts)
    print(f"error-corrected clusters (Hamming <= {CLUSTER_MAX_DIST}, >= {CLUSTER_RATIO} x abundance ratio): {len(clusters)}")
    print(f"--------------- top 10 variants (error-corrected) -----------------")
    
    for rep, members in clusters[:10]: # must cover > 80 % of the sequences 
        count = sum(c for _, c in members) # counting no. of reads in the top 10 variants
        percentage = (count/total)*100 # percentage of x sequence that exists in the sequencing data
                                         
        print(f" {count:4d} reads ({percentage:5.1f}%) {rep}")
        print(f"-> {len(members)} raw variant(s) merged")
    with open(args.valid_var_txt, "w") as f: # clustered
        for rep, members in clusters: # all the sequences availsable
            count = sum(c for _, c in members) 
            
            percentage = (count/total)*100 # percentage of x sequence that exists in the sequencing data
            f.write(f"Sequence: {rep}; percentage out of valid variants {total}: {percentage:.1f} % \n")
    print(f"-------------------------------------------------------------")
    
    with open(args.dna_var_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["read_id", "strand", "nt", "qual", "source"])
        w.writeheader()
        w.writerows(asdict(r) for r in reads)

    
    with open(args.aa_var_csv, "w", newline="") as file: # unclustered
        write = csv.DictWriter(file, fieldnames=["aa_variant", "count"])
        write.writeheader()
        for aa_seq, count in metrics.aa_counts.items():
            write.writerow({"aa_variant": aa_seq, "count":count})

    with open(args.invalid_var_discarded_at_find_anchor_csv, "w", newline="") as file_:
        write = csv.DictWriter(file_, fieldnames=["read_id", "seq", "qual", "strand", "reason"])
        write.writeheader()
        write.writerows(invalid) 
    
    with open(args.discarded_at_degen_mismatch, "w", newline="") as fi:
        write = csv.DictWriter(fi, fieldnames=["read_id", "dna_seq", "reason"])
        write.writeheader()
        write.writerows(discarded) 

if __name__ == "__main__":
    main()
 




            


    



        
        









