import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    parser = argparse.ArgumentParser()
    parser.add_argument("-start", help="start base index",required=True,type=int)
    parser.add_argument("-end", help="end base index",required=True,type=int)
    args = parser.parse_args()
    sub_seq = record.seq[args.start-1:args.end]
    print(sub_seq)
