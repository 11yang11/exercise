import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    parser = argparse.ArgumentParser()
    parser.add_argument("-start", help="start base index",required=True,type=int)
    parser.add_argument("-end", help="end base index",required=True,type=int)
    parser.add_argument("-strand", help="the output direction of sequence",choices=["+","-"])
    args = parser.parse_args()
    sub_seq = record.seq[args.start-1:args.end]
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
    else:
        print("NO parameter '-strand',please choose direction of the sequence")

