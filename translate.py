#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

#诺禾致远考试题目

import sys

def DNA_complement(nt):   #定义互补函数
    nt = nt.upper()
    nt = nt.replace('A', 't')
    nt = nt.replace('T', 'a')
    nt = nt.replace('C', 'g')
    nt = nt.replace('G', 'c')
    return nt.upper()

#构建氨基酸密码子表的字典
AAcode = {'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
'TAC':'Y','TAT':'Y','TAA':'_','TAG':'_',
'TGC':'C','TGT':'C','TGA':'_','TGG':'W'}
with open(sys.argv[1],"r") as f: 
    SeqFile = f.readlines()
    
OutFile = open(sys.argv[2],'w')

number = 0 
for seq in SeqFile:
    if seq.startswith(">"): # 提取fasta序列名称
        number = number + 1  # 记录第几条序列
        OutFile.writelines(seq) # 写入序列名称
    else:
        if number%2 == 1: # 判断奇偶数
            seq = seq[60:301] # 提取 61-300区域
            seq = ''.join([AAcode.get(seq[3*i:3*i+3],'X') for i in range(len(seq)//3)]) #将该区域翻译为氨基酸
        else:
            seq = seq.strip()
            seq = DNA_complement(seq)[::-1]  #序列反向互补
        OutFile.writelines(seq+"\n")
OutFile.close()
