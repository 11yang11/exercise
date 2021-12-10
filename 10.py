import sys
import re
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    seq = re.sub('TACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTAAAGCTAAC',
                 'TACCCTATAGCACCCCCTCTGCCCCCTCTAGAGCCCACTGTAAAGCTAAC', seq, 1)
parser = argparse.ArgumentParser()
parser.add_argument("-N",help="overlapping length between two sequences with equal length",type=int)
args = parser.parse_args()
l = len(seq)
n = args.N

def edit_distance(str1, str2):
    '''字符串str1和str2之间的最短编辑距离'''
    str1_length = len(str1)
    str2_length = len(str2)
    matrix = [[j for j in range(len(str2) + 1)] for i in range(len(str1) + 1)]
    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if str1[i-1] == str2[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = matrix[i-1][j-1]+1
    return matrix[-1][-1]

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
            diff = edit_distance(subseq_a,subseq_b)
            if diff <= 1:
                print("overlap:", subseq_a, start_a, "-", start_a + n - 1)
            start_a += 1
        start += 1
else:
    print("No parameter '-N'")
