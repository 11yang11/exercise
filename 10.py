import sys
import levenshtein
import re
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    seq = re.sub('TACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAAC',
                 'TACCCTATAGCACCCCCTCTGCCCCCTCTAGAGCCCACTGTAAAGCTAAC', seq, 1)
    h10 = (l - 51) / 2
    a = seq[:int(h10 + 50)]
    b = seq[int(h10):]
    start = 0
    while start < len(b) - 49:
        subseq_b = b[start:(start + 50)]
        start_a = 0
        while start_a < len(a) - 49:
            subseq_a = a[start_a:(start_a + 50)]
            diff = levenshtein.hamming(subseq_a, subseq_b)
            if diff <= 1:
                print("overlap:", subseq_a, start_a, "-", start_a + 49)
            start_a += 1
        start += 1
