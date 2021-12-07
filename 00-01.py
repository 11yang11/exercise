#!/usr/bin/env python
from Bio import SeqIO
for record in SeqIO.parse("human_mitochondrial.fasta", "fasta"):
    #01
    print("seq_len:", len(record.seq))
    print("gc_proportion:", (record.seq.count("C") + record.seq.count("G")) / len(record.seq))
