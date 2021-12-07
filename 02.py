import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    parser = argparse.ArgumentParser()
    parser.add_argument("-start", help="start base index", default=0,type=int)
    parser.add_argument("-end", help="end base index", default=0,type=int)
    args = parser.parse_args()
    sub_seq = record.seq[args.start:args.end]
    print(sub_seq)
