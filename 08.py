import sys
import argparse
from Bio import SeqIO
parser = argparse.ArgumentParser()
parser.add_argument("-N",help="overlapping length between two sequences with equal length",type=int)
args = parser.parse_args()
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
l = len(seq)
l = len(seq)
n = args.N
if args.N:
    if l % 2 == 0 and n % 2 == 0:
        h = (l - n) / 2
    if l % 2 != 0 and n % 2 != 0:
        h = (l - n) / 2
    else:
        n = n+1
        h = (l - n) / 2
    print(seq[:int(h + n)], '\n' * 2, seq[int(h):])
else:
    print("No parameter '-N'")
