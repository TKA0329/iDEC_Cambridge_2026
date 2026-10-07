"""
First, count length of sequence
2nd, count frequency of a certain five_prime_motif
3rd, group sequence according to length and frequency of five_prime_motif

Usage: 
    python analyzer.py input.fastq output.csv metrics.txt 
"""
from dataclasses import dataclass, field
from collections import Counter, defaultdict
import argparse
import csv
from tqdm import tqdm
import sys
from distance import l_matrix
import streamlit as st

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

def main():
    parser = argparse.ArgumentParser(
        description = "Convert a FASTQ file into a csv.")
    parser.add_argument("input_fastq", help="Path to the input .fastq file")
    parser.add_argument("output_csv", help="Path to write the output .csv file")
    parser.add_argument("metrics_txt", help="Path to write the output global metrics file")
    
    args = parser.parse_args()
    # CSV file containing read_id, seq, quality, length, and average quality of each sequence
    desired_sequence_count = 0
    # input_length = int(input("Length of desired sequence: "))
    desired_sequence = input("Desired sequence: ")
    print(len(desired_sequence))
    try:
        with open(args.output_csv, "w", newline="") as out_f:
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
                    if len(seq) == len(desired_sequence):
                        desired_sequence_count += 1 # only use pass if "if" is empty
                        # degeneracy so expected some susbtitutions
                        matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                        print(f"Sequence 1: {''.join(top)}")
                        print(f"            {''.join(match)}")
                        print(f"Sequence 2: {''.join(bot)}")
                        print(f"- means insertion or deletion; * means substitution")
                    elif len(seq) == len(desired_sequence) + 1:
                        if seq.endswith("A"):
                            print(f"Additional A at end. Suspect polyadenylation.")
                            matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                            print(f"Sequence 1: {''.join(top)}")
                            print(f"            {''.join(match)}")
                            print(f"Sequence 2: {''.join(bot)}")
                            print(f"- means insertion or deletion; * means substitution")
        
                elif seq.startswith(antisense):
                    seq=reverse_complement(seq)
                    if len(seq) == len(desired_sequence):
                        desired_sequence_count += 1
                        matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                        print(f"Sequence 1: {''.join(top)}")
                        print(f"            {''.join(match)}")
                        print(f"Sequence 2: {''.join(bot)}")
                        print(f"- means insertion or deletion; * means substitution")
                    if len(seq) == len(desired_sequence) + 1:
                        if seq.endswith("A"):
                            print(f"Additional A at end. Suspect polyadenylation.")
                            matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                            print(f"Sequence 1: {''.join(top)}")
                            print(f"            {''.join(match)}")
                            print(f"Sequence 2: {''.join(bot)}")
                            print(f"- means insertion or deletion; * means substitution")

                else: 
                    print(f"-------------------------------------------------")
                    print(f"{read_id}-> unknown sequence:")
                    print(f"{seq} \n")
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


