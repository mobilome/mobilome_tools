#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   fasta_slop.py
@Time    :   2021/12/31 21:35:20
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com
'''

# here put the import lib
import os
import re
import sys
import time
import argparse

def HowManyTime(tbegin,tend):
    """
    to calculate the time to evaluate the speed
    """
    tTotal=tend-tbegin
    tsec=tTotal%60
    ttolmin=tTotal//60
    thour=ttolmin//60
    tmin=ttolmin%60
    suretime="运行时间：%d hour, %d minutes, %.2f seconds"%(thour,tmin,tsec)
    return suretime

def write2file(name,content_list):
    with open(name,'w') as handle:
        handle.writelines(content_list)

def fasta2bed(fasta):
    bedRecords = []
    with open(fasta) as handle:
        for seq in re.findall('>.+',handle.read()):
            chrom = seq[1:].split(':')[0]
            chromStart = seq.split(':')[1].split('-')[0]
            chromEnd = seq.split('-')[1].split('(')[0].strip()
            name = '.'
            score = '.'
            strand = seq.split(')')[0].strip()[-1]
            bedRecord = '%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,chromStart,chromEnd,name,score,strand)
            bedRecords.append(bedRecord)
        return bedRecords

def bedtools_slop(bedRecords,left_len,right_len,genome_info):
    fasta_slop_bed_name = '%s-L%s-R%s.bed' % (bedRecords.split('.')[0],left_len,right_len)
    command = f'bedtools slop -s -l {left_len} -r {right_len} -g {genome_info} -i {bedRecords} > {fasta_slop_bed_name}'
    print('运行命令：%s' % command)
    os.system(command)
    return fasta_slop_bed_name

def bedtools_getfasta(fasta_input,fasta_output,fasta_slop_bed_name):
    command = f'bedtools getfasta -s -fi {fasta_input} -fo {fasta_output} -bed {fasta_slop_bed_name}'
    print('运行命令：%s' % command)
    os.system(command)


def main():
    parser = argparse.ArgumentParser( description='根据fasta文件名称提取基因组位置信息，然后前后延申并提取延申后的序列')
    parser.add_argument('-fi', '--fasta_input', type=str, required=True, help='Input FASTA file' )
    parser.add_argument('-fd', '--fasta_dataset', type=str, required=True, help='Input FASTA file' )
    parser.add_argument('-fo', '--file_output', type=str, required=True, help='Output file' )
    parser.add_argument('-l', '--left_len', type=str, required=True, help='左边延申的长度' )
    parser.add_argument('-r', '--right_len', type=str, required=True, help='右边延申的长度' )
    parser.add_argument('-g', '--genome_info', type=str, required=True, help='基因组染色体长度信息' )

    args = parser.parse_args()

    fasta_bed_name = args.fasta_input + '.bed'
    write2file(fasta_bed_name,fasta2bed(args.fasta_input))
    fasta_slop_bed_name = bedtools_slop(fasta_bed_name,args.left_len,args.right_len,args.genome_info)
    bedtools_getfasta(args.fasta_dataset,args.file_output,fasta_slop_bed_name)



if __name__ == "__main__":
    main()



