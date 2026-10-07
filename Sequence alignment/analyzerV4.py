"""
First, count length of sequence
2nd, count frequency of a certain five_prime_motif
3rd, group sequence according to length and frequency of five_prime_motif

Usage: 
    python analyzerV4.py input.fastq output.csv metrics.txt alignment.txt
"""
from dataclasses import dataclass, field
from collections import Counter, defaultdict
import argparse
import csv
from tqdm import tqdm
import sys
from distance import l_matrix



FIELDNAMES = ["read_id", "Sequence", "Quality", "Length", "Quality_average"]

@dataclass # @dataclass automatically creates the constructor for this data-holding class.
class MisprimeMetrics:
    total_reads: int=0
    length_count: Counter = field(default_factory=Counter)
    five_prime_motifs: Counter = field(default_factory=Counter)
    three_prime_junctions: defaultdict = field(default_factory=lambda: defaultdict(Counter))
    quality_sequence: dict = field(default_factory=Counter)
    quality_position: list = field(default_factory=list)
    same_3_and_5_motifs: defaultdict = field(default_factory=lambda: defaultdict(Counter))

from fastq_parser import parse_fastq

metrics = MisprimeMetrics()

quality_dict = {}
sense="AGGAAGCCAGACAAATGGTAGT"
antisense = "CCCAAGCTTCGCTCGAGCTTTA"
def calculate_stats(read_id, sequence, quality):
    length = len(sequence)
    length_of_sense= len(sense)
    length_of_antisense = len(antisense)
    motif_5 = sequence[:length_of_sense]
    motif_3 = sequence[-length_of_antisense:]

    metrics.total_reads += 1 # total number of reads 
    metrics.length_count[length] += 1 # number of sequences of each specific length 
    metrics.five_prime_motifs[motif_5] += 1 # number of sequences with a specific five prime motif
    metrics.three_prime_junctions[length][motif_3] += 1
    # number of sequences of a specific length with a specific motif_3
    metrics.same_3_and_5_motifs[motif_5][motif_3] += 1

    # converting phred scores
    convert_list = []
    for char in quality:
        convert = ord(char) -33 
        convert_list.append(convert)

    # appending converted phred score to a dict {read_id : converted list of phred socres}
    quality_dict[read_id] = convert_list
    # summing up quality of phred scores of the entire sequence
    sum_quality=sum(quality_dict[read_id])

    # putting it in a dict so {read_id: sum qaulity}
    # shud use read_id instead of sequence bcz will get overwritten as the sequence is the same
    metrics.quality_sequence[read_id] = sum_quality 

    # accummulate quality by position
    for i, score in enumerate(convert_list):
        if i >= len(metrics.quality_position):
            metrics.quality_position.append(0) # implement the check when its first created
            # prevent "list index out of range"
            # first iteration, turns [] -> [0]
            # 2nd iteration, turns [30] -> [30, 0]
            # 3rd iteration, turns [30, 32] -> [30, 32, 0]
        metrics.quality_position[i] += score # returns a final list

    quality_average= sum_quality / length
    # keeps adding score at each position i 

    return {
        "read_id":read_id,
        "Sequence": sequence,
        "Length" : length,
        "Quality": quality,
        "Quality_average": quality_average
    }
# metrics don't need to be returned bcz its already stored in global metrics object

def write_metrics(metrics, filepath):
    with open(filepath, "w") as f:
        f.write(f"Total reads: {metrics.total_reads}\n\n")
        
        f.write("Length counts: \n")
        for length, count in metrics.length_count.items():
            f.write(f"{length} bp: {count} \n")
        f.write("---------------------------------------------------")
        f.write("\n5' motifs: \n")
        f.write("\n Correct sense: 5'-AGGAAGCCAGACAAATGGTAGT-3' \n Correct antisense: 5'-CCCAAGCTTCGCTCGAGCTTTA-3' \n")
        for motif, count in metrics.five_prime_motifs.items():
            f.write(f"{motif} : {count} \n")
        f.write("---------------------------------------------------")
        f.write("\n3' junctions:\n")
        for length, motifs in metrics.three_prime_junctions.items():
            f.write(f"{length} bp: \n")
            for motif, count in motifs.items():
                f.write(f"{motif}: {count} \n")
        f.write("---------------------------------------------------")
        f.write("\n Correct sense: 5'-AGGAAGCCAGACAAATGGTAGT-3' \n Correct antisense: 5'-CCCAAGCTTCGCTCGAGCTTTA-3' \n")
        f.write("\n5'/3' motif combinations:\n")
        for five, three in metrics.same_3_and_5_motifs.items():
            f.write(f"5'-> {five}: \n")
            for motif, count in three.items():
                f.write(f"      3'->{motif}: {count} \n")
        f.write("---------------------------------------------------")
        f.write("\nQuality sums by position:\n")
        for position, total in enumerate(metrics.quality_position):
            f.write(f"Position {position}: {total}\n")


# reverse_complement
def reverse_complement(sequence):
    complement = {"A":"T",
                  "T":"A",
                  "C":"G",
                  "G":"C",
                  "N": "N"}
    return "".join(complement[base] for base in reversed(sequence))

def complement(sequence):
    complement = {"A":"T",
                "T":"A",
                "C":"G",
                "G":"C",
                "N": "N"}
    return "".join(complement[base] for base in sequence)

def print_clustal_style(aln_f, top, match, bot, width=60, name1="Sequence 1", name2="Sequence 2"):
    """
    Takes the three aligned strings + 2 optional settings: width
    (how many alignment columns to show per line-block, defaulting to 60) + the two display names for the sequences
    """
    # don't do with oepn file here cz will get rerun each time in the loop in main()
    # top, match, bot are all the same length so len(top) tells you the total alignment width
    aln_len = len(top)

    # counters for actual pos in each o.g. sequence
    # how far am I into the actual biological sequence, ignoring the bookkeeping dashes the alignment added?
    pos1=0 # real (non-gap) base count consumed so far in top 
    pos2=0 # real (non-gap) base count consumed so far in bot 

    label_width= max(len(name1), len(name2)) # keeps alignment bars starting at same col

    # automatically produces however many chunks needed 
    for start in range(0, aln_len, width):
        # last chunk could overshoot, min caps it at aln_len so the last chunk comes out shorter
        end = min(start+width, aln_len) 

        # all top, match, bot columns are aligned
        top_chunk = top[start:end]
        match_chunk = match[start:end]
        bot_chunk=bot[start:end]

        # advance position counters, ignoring gap chars, counts how many top chars are not a gap char("-") 
        pos1 += sum(1 for c in top_chunk if c != "-")
        pos2 += sum(1 for c in bot_chunk if c != "-")

        #aln_f.write(f"---------------------------------------------------------- \n")
        aln_f.write(f"{name1:<{label_width}}  {top_chunk}  {pos1} \n")
        aln_f.write(f"{'':<{label_width}}  {match_chunk} \n")
        aln_f.write(f"{name2:<{label_width}}  {bot_chunk}  {pos2} \n")
        aln_f.write(f"\n") # need to receive 1 pos argument
        # aln_f.write(f"----------------------------------------------------------- \n")

def main():
    parser = argparse.ArgumentParser(
        description = "Convert a FASTQ file into a csv.")
    parser.add_argument("input_fastq", help="Path to the input .fastq file")
    parser.add_argument("output_csv", help="Path to write the output .csv file")
    parser.add_argument("metrics_txt", help="Path to write the output global metrics file")
    parser.add_argument("alignment_txt", help="Path to write the output alignment file")
    args = parser.parse_args()
    # CSV file containing read_id, seq, quality, length, and average quality of each sequence
    desired_sequence_count = 0
    # input_length = int(input("Length of desired sequence: "))
    desired_sequence = input("Desired sequence: ")
    print(f"Length of seuquence: {len(desired_sequence)}")
    try:
        with open(args.output_csv, "w", newline="") as out_f, \
            open(args.alignment_txt, "w") as aln_f:
            writer = csv.DictWriter(out_f, fieldnames=FIELDNAMES)
            writer.writeheader()
            for read_id, seq, qual in tqdm(
                parse_fastq(args.input_fastq),
                desc = "Processing reads",
                unit = "read", 
            ):
                writer.writerow(calculate_stats(read_id, seq, qual))
                # if len(seq) == input_length:
                # sense primer       → keep
                # antisense primer   → reverse complement
                # neither            → flag as unknown
                
                if seq.startswith(sense):
                    aln_f.write(f"------------------------------------------------- \n")
                    if len(seq) == len(desired_sequence) + 1:
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r}; ENDS WITH: {seq[-25:]!r} \n")
                        if seq.endswith("A"):
                            aln_f.write(f"{read_id} -> Additional A at end. Possible polyadenylation. \n")
                            matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                            print_clustal_style(aln_f, "".join(top), "".join(match), "".join(bot),
                    name1="original", name2="pcr product")
                        else:
                            aln_f.write(f"{read_id} -> Correct 5' motif, +1 length, but not polyA (possible insertion): \n")
                            aln_f.write(f"{seq}\n")
                    elif len(seq) == len(desired_sequence) and seq.endswith(reverse_complement(antisense)):
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r} \n")
                        aln_f.write(f"{read_id} -> correct sequence? Check for indels. \n")
                        desired_sequence_count += 1 # only use pass if "if" is empty
                        # degeneracy so expected some susbtitutions
                        matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                        print_clustal_style(aln_f, "".join(top), "".join(match), "".join(bot),
                    name1="original", name2="pcr product")
                    

                    else:
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r};  ENDS WITH: {seq[-25:]!r} \n")
                        aln_f.write(f"{read_id} -> Correct 5' motif but incorrect length and incorrect 3' motif: \n")
                        aln_f.write(f"{seq}\n")

                elif seq.startswith(antisense):
                    aln_f.write(f"------------------------------------------------- \n")
                    aln_f.write(f"NOTE: This sequence was reverse complemented before analysis \n")
                    seq=reverse_complement(seq)

                    if len(seq) == len(desired_sequence) + 1:
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r}; ENDS WITH: {seq[-25:]!r} \n")
                        if seq.endswith("A"):
                            aln_f.write(f"{read_id} -> Additional A at end. Possible polyadenylation. \n")
                            matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                            print_clustal_style(aln_f, "".join(top), "".join(match), "".join(bot),
                    name1="original", name2="pcr product")
                        else:
                            aln_f.write(f"{read_id} -> Correct 5' motif, +1 length, but not polyA (possible insertion): \n")
                            aln_f.write(f"{seq}\n")

                    elif len(seq) == len(desired_sequence) and seq.startswith(sense):
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r}; ENDS WITH: {seq[-25:]!r} \n")
                        aln_f.write(f"{read_id} -> correct sequence? Check for indels. \n")
                        desired_sequence_count += 1
                        matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                        print_clustal_style(aln_f, "".join(top), "".join(match), "".join(bot),
                    name1="original", name2="pcr product")
                    
                    else:
                        aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r}; ENDS WITH: {seq[-25:]!r} \n")
                        aln_f.write(f"{read_id} -> Correct 5' motif but incorrect length and incorrect 3' motif: \n")
                        aln_f.write(f"{seq}\n")

                else: 
                    aln_f.write(f"------------------------------------------------- \n")
                    aln_f.write(f"RAW SEQ STARTS WITH: {seq[:25]!r}; ENDS WITH: {seq[-25:]!r} \n")
                    aln_f.write(f"{read_id} -> unknown sequence: \n")
                    aln_f.write(f"Does not start with the correct motifs")
                    aln_f.write(f"{seq} \n")
                    matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                    print_clustal_style(aln_f, "".join(top), "".join(match), "".join(bot),
                                        name1="original", name2="pcr product")
                    # aln_f.write(f"{'-'*50} \n")
                    continue # if unknown orientation, stop further analysis 
                
            print(f"Number of desired sequences in sequencing data:  {desired_sequence_count}")

                    
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    # Txt file containing the metrics above (e.g. how many of them have a certain length)
    try:
        write_metrics(metrics, args.metrics_txt)
    except ValueError as err:
        print(f"Error: {err}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Done. Wrote stats for {args.input_fastq} to {args.output_csv} and {args.metrics_txt}.")

if __name__=="__main__":
    main()