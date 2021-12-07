import sys
import argparse
import re
import levenshtein
#00
sys.path.append("d:\lib\site-packages")
from Bio import SeqIO
for record in SeqIO.parse("d:\BGI\human_mitochondrial.fasta", "fasta"):
    #01
    print("seq_len:", len(record.seq))
    print("gc_proportion:", (record.seq.count("C") + record.seq.count("G")) / len(record.seq))
    parser = argparse.ArgumentParser()
    parser.add_argument("-start", help="start base index", default=0,type=int)
    parser.add_argument("-end", help="end base index", default=0,type=int)
    parser.add_argument("-g", help="gc_content of subseqence",action="store_true")
    parser.add_argument("-strand", help="the output direction of sequence",choices=["+","-"])
    parser.add_argument("-window",help="the width of window",type=int)
    parser.add_argument("-step",help="the steps every time",type=int)
    parser.add_argument("-n",help="the number of overlapping bases",type=int)
    args = parser.parse_args()
    #02
    sub_seq = record.seq[args.start:args.end]
    print(sub_seq)
    #03
    if args.g:
        print((sub_seq.count("C")+sub_seq.count("G"))/len(sub_seq))
    #04
    if args.strand == "+":
        print(sub_seq)
    elif args.strand == "-":
    #方法1（dictionary）
        complement = {"A":"T","T":"A","C":"G","G":"C"}
        sub_seq_list = list(sub_seq)
        complement_list = [complement[base] for base in sub_seq_list]
        #遍历子序列列表，将列表中的碱基赋值给字典键
        recerse_com = ''.join(complement_list[::-1])
        #把列表中元素转化成字符串
        print(recerse_com)
    #方法2（sequence）
        sub_seq = sub_seq.replace('A','t')
        sub_seq = sub_seq.replace('T','a')
        sub_seq = sub_seq.replace('C','g')
        sub_seq = sub_seq.replace('G','c')
        #如果先将字符串A转换成T，再将T转换成A的时候就会出问题，此时原始的T和由A转换的T都会被转换掉。
        rcomplement_sub_seq = sub_seq.upper()
        print(rcomplement_sub_seq[::-1])
    #方法3（translate)
        sub_seq = ''.join(sub_seq)
        #translate函数只能在字符串中使用，必须把上面的sub_seq'Bio.Seq.Seq‘格式转化为字符串形式
        transtable = sub_seq.maketrans('ATCG','TAGC')
        complement_seq = sub_seq.translate(transtable)
        print(complement_seq[::-1])
    #05
    window = args.window
    step = args.step
    for index in range(len(record.seq)):
        if index % step == 0:
            print(record.id, "_s", index, "e", index + window, ":", record.seq[index:index + window], sep='')
    #06
    seq = ''.join(record.seq)
    transtable = seq.maketrans('ATCG', 'TAGC')
    complement_seq = seq.translate(transtable)
    re_com_seq = complement_seq[::-1]
    # translate方法得到反向互补序列
    forward_seq = re.search('catcaactga', seq, re.I).span()
    reverse_seq = re.search('catcaactga', re_com_seq, re.I).span()
    print('"catcaactga"forward_seq:', forward_seq, '\n"catcaactga"reverse_seq:', reverse_seq)
    # 07
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
    # 08
    l = len(seq)
    if l % 2 == 0:
        n = 10
    else:
        n = 11
    h = (l - n) / 2
    print(seq[:int(h + n)], '\n' * 2, seq[int(h):])
    # 09
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
    # 10
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

