"""
First, count length of sequence
2nd, count frequency of a certain five_prime_motif
3rd, group sequence according to length and frequency of five_prime_motif

Usage: 
    python analyzer.py input.fastq output.csv metrics.txt 
"""
import tempfile
from dataclasses import dataclass, field
from collections import Counter, defaultdict
import argparse
import csv
from tqdm import tqdm
import sys
from distance import l_matrix
import streamlit as st

st.title("Sequencing data analyzer.")
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
import io
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

def write_metrics(metrics):
    buf = io.StringIO()
    buf.write(f"Total reads: {metrics.total_reads}\n\n")
    
    buf.write("Length counts: \n")
    for length, count in metrics.length_count.items():
        buf.write(f"{length} bp: {count} \n")
    buf.write("---------------------------------------------------")
    buf.write("\n5' motifs: \n")
    buf.write("\n Correct sense: 5'-AGGAAGCCAGACAAATGGTAGT-3' \n Correct antisense: 5'-CCCAAGCTTCGCTCGAGCTTTA-3' \n")
    for motif, count in metrics.five_prime_motifs.items():
        buf.write(f"{motif} : {count} \n")
    buf.write("---------------------------------------------------")
    buf.write("\n3' junctions:\n")
    for length, motifs in metrics.three_prime_junctions.items():
        buf.write(f"{length} bp: \n")
        for motif, count in motifs.items():
            buf.write(f"{motif}: {count} \n")
    buf.write("---------------------------------------------------")
    buf.write("\n Correct sense: 5'-AGGAAGCCAGACAAATGGTAGT-3' \n Correct antisense: 5'-CCCAAGCTTCGCTCGAGCTTTA-3' \n")
    buf.write("\n5'/3' motif combinations:\n")
    for five, three in metrics.same_3_and_5_motifs.items():
        buf.write(f"5'-> {five}: \n")
        for motif, count in three.items():
            buf.write(f"      3'->{motif}: {count} \n")
    buf.write("---------------------------------------------------")
    buf.write("\nQuality sums by position:\n")
    for position, total in enumerate(metrics.quality_position):
        buf.write(f"Position {position}: {total}\n")
    

# reverse_complement
def reverse_complement(sequence):
    complement = {"A":"T",
                  "T":"A",
                  "C":"G",
                  "G":"C",
                  "N": "N"}
    return "".join(complement[base] for base in reversed(sequence))

def print_clustal_style(buf, top, match, bot, width=60, name1="Sequence 1", name2="Sequence 2"):
    """
    Takes the three aligned strings + 2 optional settings: width
    (how many alignment columns to show per line-block, defaulting to 60) + the two display names for the sequences
    """
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

        buf.write(f"{name1:<{label_width}}  {top_chunk}  {pos1}")
        buf.write(f"{'':<{label_width}}  {match_chunk}")
        buf.write(f"{name2:<{label_width}}  {bot_chunk}  {pos2}")
    buf.write("\n")

def main():
    uploaded_fastq = st.file_uploader("Upload FASTQ File: ", type = ["fastq", "fq"])

    # CSV file containing read_id, seq, quality, length, and average quality of each sequence
    desired_sequence_count = 0
    # input_length = int(input("Length of desired sequence: "))
    desired_sequence = st.text_input("Desired sequence: ")
    # print(len(desired_sequence))
    if uploaded_fastq and desired_sequence and st.button("Run analysis"):
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".fastq", delete=False) as tmp:
            tmp.write(uploaded_fastq.getvalue())
            tmp_path = tmp.name
        csv_buf = io.StringIO()
        aln_buf = io.StringIO()

        writer = csv.DictWriter(csv_buf, fieldnames=FIELDNAMES)
        writer.writeheader()
    
        for read_id, seq, qual in parse_fastq(tmp_path):
            writer.writerow(calculate_stats(read_id, seq, qual))
            # if len(seq) == input_length:
            # sense primer       → keep
            # antisense primer   → reverse complement
            # neither            → flag as unknown
            st.subheader("Alignment")
            if seq.startswith(sense):
                if len(seq) == len(desired_sequence):
                    desired_sequence_count += 1 # only use pass if "if" is empty
                    # degeneracy so expected some susbtitutions
                    matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                    print_clustal_style(aln_buf, "".join(top), "".join(match), "".join(bot),
                name1="original", name2="pcr product")
                elif len(seq) == len(desired_sequence) + 1:
                    if seq.endswith("A"):
                        #st.write(f"Additional A at end. Suspect polyadenylation.")
                        matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                        print_clustal_style(aln_buf, "".join(top), "".join(match), "".join(bot),
                name1="original", name2="pcr product")
            elif seq.startswith(antisense):
                seq=reverse_complement(seq)
                if len(seq) == len(desired_sequence):
                    desired_sequence_count += 1
                    matrix, operations, top, match, bot = l_matrix(desired_sequence, seq)
                    print_clustal_style(aln_buf, "".join(top), "".join(match), "".join(bot),
                name1="original", name2="pcr product")
                if len(seq) == len(desired_sequence) + 1:
                    if seq.endswith("A"):
                        # st.write(f"Additional A at end. Suspect polyadenylation.")
                        matrix, operations, top, match, bot = l_matrix(desired_sequence[:len(seq)-1], seq)
                        print_clustal_style(aln_buf, "".join(top), "".join(match), "".join(bot),
                name1="original", name2="pcr product")

        st.write(f"Number of desired sequences in sequencing data:  {desired_sequence_count}")

    # Txt file containing the metrics above (e.g. how many of them have a certain length)
        metrics_str = write_metrics(metrics)
        
        with st.expander("View alignments"):
            st.text(aln_buf.getvalue())
        st.download_button("Download CSV", data=csv_buf.getvalue(), file_name="output.csv", mime="text/csv")
        st.download_button("Download alignments", data=aln_buf.getvalue(), file_name="alignments.txt", mime="text/plain")
        st.download_button("Download metrics", data=metrics_str, file_name="metrics.txt", mime="text/plain")

if __name__=="__main__":
    main()