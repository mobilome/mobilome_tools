#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#诺禾致远考试题目

import sys,re

with open(sys.argv[1],"r") as f: 
    SeqFile = f.readlines()
Sign_Chr = sys.argv[2]
#print("{0}\t{1}\t{2}\t{3}\t{4}\n" .format("Seq_ID","Seq_Length(bp)","GC%","Sign_num","Sign_pos"))

for seq in SeqFile:
    if seq.startswith(">"): # 提取fasta序列名称
        seq  = seq.strip()
        Seq_ID = seq[1:len(seq)]
        #print(Seq_ID,end='\t')
    else:
        Seq_Length = len(seq)
        seq_GC = len(re.findall('[GCgc]', seq))
        GC_precent = round(seq_GC/Seq_Length*100,2)
        Sign_num = seq.count(Sign_Chr)
        Sign_pos = []
        for a in re.finditer(Sign_Chr, seq):
            print(a.span())
            #Sign_pos.append(a.span())
        #print("{}\t{}\t{}\t{}\n" .format(Seq_Length,GC_precent,Sign_num,Sign_pos))

