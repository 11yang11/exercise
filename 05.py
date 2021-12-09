import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    parser = argparse.ArgumentParser()
    parser.add_argument("-window",help="the width of window",required=True,type=int)
    parser.add_argument("-step",help="the steps every time",required=True,type=int)
    args = parser.parse_args()
    window = args.window
    step = args.step
    for index in range(len(record.seq)):
        if len(record.seq)-index >= window and index % step == 0:
            print(record.id, "_s", index, "e", index + window, ":", record.seq[index:index + window], sep='')
