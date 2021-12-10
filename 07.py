import sys
import re
from Bio import SeqIO

def edit_distance(str1, str2):
'''字符串str1和str2之间的最短编辑距离'''
    str1_length = len(str1)
    str2_length = len(str2)
    matrix = [[j for j in range(len(str2) + 1)] for i in range(len(str1) + 1)]
    #初始化一个二维数组
    for i in range(1,len(str1)+1):
        for j in range(1,len(str2)+1):
            if str1[i-1] == str2[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = matrix[i-1][j-1]+1
    return matrix[-1][-1]
    #返回最后一个值
    
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)    
    for index in range(len(seq)):
        if len(seq) - index > 10:
            subseq = seq[index:index + 10]
            subseq = subseq.lower()
        difference = edit_distance(subseq, 'catcaactta')
        if difference <= 1:
            print(re.search(subseq, seq, re.I), '\n')
