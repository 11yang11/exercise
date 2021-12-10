import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
parser = argparse.ArgumentParser()
parser.add_argument("-N",help="overlapping length between two sequences with equal length",type=int)
args = parser.parse_args()
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
    a = seq[:int(h + n)]
    b = seq[int(h):]
    start = 0
    while start < len(b)- n + 1:
        subseq_b = b[start:(start + n)]
        start_a = 0
        while start_a < len(a) - n + 1:
            subseq_a = a[start_a:(start_a + n)]
            if subseq_a == subseq_b:
                print("overlap:", subseq_a, start_a, "-", start_a + n - 1)
            start_a += 1
        start += 1
else:
    print("No parameter '-N'")
