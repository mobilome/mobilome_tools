#!/usr/bin/env python3

# -*- encoding: utf-8 -*-

"""
@Author  :   Naisu
 
@Contact :   3298990@qq.com
 
This script can be used to transform RepeatMasker out file into a BED file.
"""

import argparse


def main():

    ## -i and -o for input/output.
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
	parser.add_argument('-r', '--repeat_class', type=str, required=True, help='pick up special repeat_class' )
    args = parser.parse_args()

    with open(args.input_file) as f:
        file = f.readlines()
    opf = open(args.output_file,'w')
    
    for line in file:
        col = line.split()
        if len(col) > 14:
            SW_score = col[0]
            perc_div = col[1]
            perc_del = col[2]
            perc_ins = col[3]
            query_sequence = col[4]
            piq_begin = col[5]
            piq_end = col[6]
            piq_left = col[7]
            strand = col[8]
            match_repeat = col[9]
            repeat_family = col[10]
            pir_begin = col[11]
            pir_end = col[12]
            pir_left = col[13]
            ID = col[14]
            if repeat_family == args.repeat_class:
                if strand == "C":
                    strand = "-"
                opf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n" .format(chrom,chromStart,chromEnd,'Default','0',strand))

    opf.close()

if __name__ == '__main__':
    main()
