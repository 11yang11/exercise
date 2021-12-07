import sys
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    l = len(seq)
    if l % 2 == 0:
        n = 10
    else:
        n = 11
    h = (l - n) / 2
    print(seq[:int(h + n)], '\n' * 2, seq[int(h):])
