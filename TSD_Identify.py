#!/usr/bin/env python3

# -*- encoding: utf-8 -*-

"""
@Author  :   Naisu
 
@Contact :   3298990@qq.com
 
This script can be used to identify TSD in SINE.
"""

import argparse

def Check_TSD(Five_Flank,Three_Flank,TSD_Length): #定义三个参数：5’ 序列，3’序列，TSD的长度
    TSD = []
    for index in range(len(Five_Flank)):
        Repeat_seq = Five_Flank[index:index+TSD_Length]
        if (len(Repeat_seq)) == TSD_Length and Repeat_seq in Three_Flank:
            TSD.append(Repeat_seq)
    return TSD


def main():
    parser = argparse.ArgumentParser( description='Fasta coverter to Bed format')

    ## -i and -o for input/output.
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    parser.add_argument('-f', '--Five_Flank_Length', type=int, default=100, help='Search for the length of TSD in five flank, default=100' )
    parser.add_argument('-t', '--Three_Flank_Length', type=int, default=100, help='Search for the length of TSD in three flank, default=100' )
    parser.add_argument('-l', '--Flank_Length', type=int, default=25, help='The Max length for the TSD ,default=25' )
    args = parser.parse_args()

    with open(args.input_file) as f:
        Seqs = f.readlines()
    opf = open(args.output_file,'w')
    
    for i in Seqs:
        if i.startswith(">"):
            opf.write(i)
        else:
            Five_Flank = i[:args.Five_Flank_Length]
            Three_Flank = i[-(args.Three_Flank_Length +1 ):-1]   
            for i in range(args.Flank_Length):
                TSD_Number = Check_TSD(Five_Flank,Three_Flank,i)
                if len(TSD_Number) == 1 or len(TSD_Number) == 0:
                    opf.write("".join(TSD_Number)+"\n")
                    break
                
    opf.close()

if __name__ == '__main__':
    main()
