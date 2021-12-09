import sys
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    h2 = (l - 51) / 2
    a = seq[:int(h2 + 50)]
    b = seq[int(h2):]
    start = 0
    while start < len(b) - 49:
        subseq_b = b[start:(start + 50)]
        start_a = 0
        while start_a < len(a) - 49:
            subseq_a = a[start_a:(start_a + 50)]
            if subseq_a == subseq_b:
                print("overlap:", subseq_a, start_a, "-", start_a + 49)
            start_a += 1
        start += 1