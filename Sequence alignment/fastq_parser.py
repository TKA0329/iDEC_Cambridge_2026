import argparse 
import csv
# FIELDNAMES = ["read_id", "seq", "quality"]

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
