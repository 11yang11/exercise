import sys
import re
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    seq = seq.upper()
    transtable = seq.maketrans('ATCG', 'TAGC')
    complement_seq = seq.translate(transtable)
    re_com_seq = complement_seq[::-1]
    # translate方法得到反向互补序列
    forward_seq = re.search('catcaactga', seq, re.I).span()
    reverse_seq = re.search('catcaactga', re_com_seq, re.I).span()
    print('"catcaactga"forward_seq:', forward_seq, '\n"catcaactga"reverse_seq:', reverse_seq)
