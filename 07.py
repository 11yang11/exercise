import sys
import re
import levenshtein
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    seq = ''.join(record.seq)
    #原始方法(正则表达式）
    pattern = re.compile('''.atcaactta|c.tcaactta|ca.caactta|cat.aactta|catc.actta|catca.ctta|
                 catcaa.tta|catcaac.ta|catcaact.a|catcaactt.''', re.I)
    positions = pattern.finditer(seq)
    for position in positions:
        print(position)
    # 滑动窗口
    for index in range(len(seq)):
        if len(seq) - index > 10:
            subseq = seq[index:index + 10]
            subseq = subseq.lower()
        difference = levenshtein.hamming(subseq, 'catcaactta')
        if difference <= 1:
            print(re.search(subseq, seq, re.I), '\n')
