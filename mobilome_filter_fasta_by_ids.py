#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   mobilome_filter_fasta_by_ids.py
@Time    :   2021/10/23 10:54:56
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com

此脚本文件用来从多序列文件中根据序列号或者id提取对应的子集
输入文件为fasta文件和一个包含需要提取的序列的id的单列文件
	- a text file containing a list of IDs
	- a single ID or comma separated list of IDs on the command line
'''

# here put the import lib
import os
import sys
import time
import argparse
import pandas as pd

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
    
def create_index(database_file):
    t1 = time.time()
    print('索引文件不存在,创建索引文件')
    #通过grep查询数据库对应的所有id及其对应的行号
    os.system(f"grep -n '>' '{database_file}' | cut -d ' ' -f 1 > '{database_file}_row_index_and_fasta_name'")
    #将提取的行号和序列名称分开
    os.system(f"cut -d ':' -f 1 '{database_file}_row_index_and_fasta_name' > '{database_file}_row_index'")
    os.system(f"cut -d '>' -f 2 '{database_file}_row_index_and_fasta_name' > '{database_file}_fasta_name'")

    #将行号删除第一行，并在最后一行加上数据库的总行数
    os.system(f"sed '1d' '{database_file}_row_index' > '{database_file}_row_index_shift'")
    os.system(f"wc -l '{database_file}' | cut -d ' ' -f 1  >> '{database_file}_row_index_shift'")

    #构建出数据库每条序列所对应的行号
    os.system(f"paste -d, '{database_file}_fasta_name' '{database_file}_row_index' '{database_file}_row_index_shift' > '{database_file}.index'")
    os.system(f"rm '{database_file}_row_index_and_fasta_name'")
    os.system(f"rm '{database_file}_row_index'")
    os.system(f"rm '{database_file}_fasta_name'")
    os.system(f"rm '{database_file}_row_index_shift'")
    print("----------完成----------")
    t2 = time.time()
    print("创建索引 %s"%HowManyTime(t1,t2))
    
def extract_fasta(database_file,filter_id_file,output_file):
    t3 = time.time()
    #通过pandas的merge函数，找出需要提取的fasta及其所在行
    print(f"开始读取索引文件")
    database = pd.read_csv(f'{database_file}.index',header=None)
    print("----------完成----------")
    print(f"开始读取待提取列表文件")
    filter_id = pd.read_csv(filter_id_file,header=None,delim_whitespace=True)
    #print(filter_id)
    print("----------完成----------")
    print("生成待提取列表索引")
    filter_id_list = filter_id[0].values.tolist()
    merge = database[database.iloc[:,0].isin(filter_id_list)]
    with open(f"{database_file}_filter_fasta_index.csv",'w') as filter_index:
        for tup in zip(merge[1],merge[2]):
            filter_index.writelines([str(i)+'\n' for i in range(int(tup[0]),int(tup[1]))])
    #merge.to_csv(f"{database_file}_filter_fasta_index.csv",index=False, header=None)
    print("----------完成----------")
    print("开始提取文件")
    awk_command = r"awk 'NR==FNR{d[$1]; next}FNR in d'" + f" '{database_file}_filter_fasta_index.csv' '{database_file}' > '{output_file}'"
    os.system(awk_command)
    os.system(f"rm '{database_file}_filter_fasta_index.csv'")
    print("----------完成----------")
    print(f"结果保存在:{output_file}")
    t4 = time.time()
    print("提取序列 %s"%HowManyTime(t3,t4))

def main():
    parser = argparse.ArgumentParser( description='从多序列文件中,根据序列号（id/accession）提取序列子集')
    parser.add_argument('-f', '--fasta_file', type=str, required=True, help='multi-FASTA文件' )
    parser.add_argument('-i', '--id_list', type=str, required=True, help='需要提取的序列id' )
    parser.add_argument('-o', '--output_file', type=str, required=True, help='输出文件名' )
    args = parser.parse_args()

    database_file = args.fasta_file #提取的目标数据库
    filter_id_file = args.id_list #需要提取序列的id
    
    if os.path.exists(f'{database_file}.index'):
        print('索引文件已存在,开始提取序列')
        extract_fasta(database_file,filter_id_file,args.output_file)
    else:
        
        create_index(database_file)
        extract_fasta(database_file,filter_id_file,args.output_file)


if __name__ == '__main__':
    main()
