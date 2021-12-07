import sys
import argparse
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="gc_content of subseqence",action="store_true")
    args = parser.parse_args()
    if args.g:
        print((sub_seq.count("C")+sub_seq.count("G"))/len(sub_seq))
