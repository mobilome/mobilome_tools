#!/usr/bin/env python3

# -*- encoding: utf-8 -*-

"""
@Author  :   Naisu
 
@Contact :   3298990@qq.com
 
This script can be used to transform a Fasta file into a BED file.

"""

import argparse


def main():
    parser = argparse.ArgumentParser( description='Fasta coverter to Bed format')

    ## -i and -o for input/output.
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to an output file to be created' )
    args = parser.parse_args()

    if args.output_file is None:
        opf = sys.stdout
    else:
        opf = open ( args.output_file, 'w')

    for line in open(args.input_file):
        if line.startswith(">"):
            chrom = line.split(">")[1].split(":")[0]
            chromStart = line.split(">")[1].split(":")[1].split("-")[0]
            chromEnd = line.split(">")[1].split(":")[1].split("-")[1].split("(")[0]
            strand = line.split(">")[1].split(":")[1].split("-")[1].split("(")[1].split(")")[0]
            if strand != '+':
                strand = '-'
        opf.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n" .format(chrom,chromStart,chromEnd,'Default','0',strand))

    opf.close()

if __name__ == '__main__':
    main()
