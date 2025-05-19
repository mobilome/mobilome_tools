#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   mobilome_tn_mining_stats.py
@Time    :   2025/04/30 10:42:27
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com
'''

# here put the import lib
import sys
import argparse
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO
import pandas as pd
from typing import List, Dict, Tuple, Optional

def calculate_statistics(lengths: List[int]) -> Tuple[int, int, int]:
    """Calculate min, max, and average of lengths."""
    if not lengths:
        return 0, 0, 0
    return min(lengths), max(lengths), round(sum(lengths) / len(lengths))

def run_cmd(cmd: str) -> None:
    try:
        subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,shell=True)
    except Exception as e:
        print(e)
        sys.exit(1)

def analyze_fasta_files(fasta_folder: Path, min_copy: Optional[int] = None, max_copy:Optional[int] = None):
    fasta_path = Path(fasta_folder)
    results = []
    
    # 遍历所有.fasta文件
    for fasta_file in fasta_path.glob('*.fasta'):
        file_name = fasta_file.name
        lengths = []

        try:
            # 使用BioPython解析FASTA文件
            for record in SeqIO.parse(fasta_file, "fasta"):
                lengths.append(len(record.seq))
            
            if lengths:  # 如果文件中有序列
                file_stats = {
                    "filename": file_name,
                    "num_sequences": len(lengths),
                    "min_length": min(lengths),
                    "max_length": max(lengths),
                    "avg_length": round(sum(lengths)/len(lengths), 2)
                }
            else:
                file_stats = {
                    "filename": file_name,
                    "num_sequences": 0,
                    "min_length": None,
                    "max_length": None,
                    "avg_length": None
                }
            results.append(file_stats)
            
        except Exception as e:
            print(f"Error processing file {file_name}: {e}")
            results.append({
                "filename": file_name,
                "num_sequences": None,
                "min_length": None,
                "max_length": None,
                "avg_length": None,
                "error": str(e)
            })
    
    # 转换为DataFrame便于查看和分析
    if results:
        df = pd.DataFrame(results)

        # 根据min_copy和max_copy过滤num_sequences
        if min_copy is not None:
            df = df[df['num_sequences'] >= min_copy]
        if max_copy is not None:
            df = df[df['num_sequences'] <= max_copy]

        df.sort_values(by='num_sequences', ascending=True, inplace=True, na_position='last')
        df = df.reset_index(drop=True)
        return df
    return pd.DataFrame()

def analyze_logs_files(logs_folder, output_file=None):
    logs_path = Path(logs_folder)
    
    # Initialize counters
    counters = {
        'total_files': 0,
        'empty_tblastn_files': 0,
        'success_files': 0
    }
    
    # Analyze all .log files
    for log_file in logs_path.glob('*.log'):
        counters['total_files'] += 1
        
        try:
            content = log_file.read_text(encoding='utf-8')
            
            # 检查是否包含特定字符串
            if "tblastn results are empty" in content:
                counters['empty_tblastn_files'] += 1
            #if "All steps completed successfully" in content:
            elif "successfully" in content:
                counters['success_files'] += 1
            else:
                print(f"\033[91m\u2717 Pay Attention!\033[0m job execution unsuccessful found in {log_file}")
        except Exception as e:
            print(f"Error reading file {log_file.name}: {e}")

    # Calculate others count
    others = counters['total_files'] - counters['empty_tblastn_files'] - counters['success_files']

    # Prepare output lines
    output_lines = [
        f"Total logs files:\t{counters['total_files']}",
        f"Empty Results:\t\t{counters['empty_tblastn_files']}",
        f"Successful Results:\t{counters['success_files']}",
        f"Others:\t\t\t{others}"
    ]

    # Output to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            print('\n'.join(output_lines), file=f)
    
    # Always print to console
    print('\n'.join(output_lines))

def run_ORFfinder(fasta_folder, results_dir: Path, args):

    fasta_df = analyze_fasta_files(fasta_folder,args.min_copy,args.max_copy)
    if fasta_df.empty:
        print("没有符合条件的结果")
        return None

    # 创建ORFfinder输出文件夹
    orf_outdir = Path(results_dir, 'ORFfinder_out')
    if orf_outdir.is_dir():
        shutil.rmtree(orf_outdir)
    orf_outdir.mkdir(exist_ok=False, parents=False)

    # 运行ORFfinder
    orfml = "orf_gr_" + str(args.min_orf) + "nt"
    orf_stats = []
    # 遍历所有.fasta文件

    for fasta_file in fasta_folder.glob('*.fasta'):
        outfile = Path(orf_outdir,fasta_file.name.replace('.fasta', '.orf.fasta'))
        rename_fasta_file = Path(orf_outdir,fasta_file.name.replace('.fasta', '.rename.fasta'))
        renameseqs = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.id, record.description = record.id.split(")", 1)
            record.id = record.id + ')'
            record.description = record.description[1:]
            if len(record.id) > 50:
                print(f'序列名称太长')
                sys.exit(1)
            renameseqs.append(record)
        SeqIO.write(renameseqs, rename_fasta_file, "fasta")
        cmd = f"ORFfinder -in {rename_fasta_file} -out {outfile} -ml {args.min_orf} -outfmt 0"
        run_cmd(cmd)
        if Path(outfile).stat().st_size == 0:
            outfile.unlink()
            continue 
        seqnames = {record.id.split("_")[1].split("(")[0] for record in SeqIO.parse(outfile, "fasta")}
        # # Find matching sequences in original file
        # sequences = []
        # for record in SeqIO.parse(fasta_file, "fasta"):
        #     if any(seqname in record.id for seqname in seqnames):
        #         sequences.append(record)
        # if sequences:
        #     filter_seq = Path(orf_outdir,fasta_file.name)
        #     SeqIO.write(sequences, filter_seq, "fasta")
        #     orf_stats.append({"filename":fasta_file.name, orfml:len(seqnames)})
        orf_stats.append({"filename":fasta_file.name, orfml:len(seqnames)})
    merge_results = pd.merge(fasta_df, pd.DataFrame(orf_stats), on='filename')
    merge_results['percent'] = (merge_results[orfml] / merge_results['num_sequences']).round(2)
    if args.min_perc:
        # 创建filter输出文件夹
        filter_seq_outdir = Path(results_dir, 'filter_seq_out')
        if filter_seq_outdir.is_dir():
            shutil.rmtree(filter_seq_outdir)
        filter_seq_outdir.mkdir(exist_ok=False, parents=False)
        filtered_df = merge_results[merge_results['percent'] >= args.min_perc]
        if filtered_df.empty:
            print("没有符合条件的结果")
            return None
        for filename in filtered_df['filename']:
            src_path = Path(fasta_folder,filename)
            dst_path = Path(filter_seq_outdir,filename)
            shutil.copy(src_path, dst_path)
        filtered_df.to_csv(Path(results_dir, 'statistics.csv'), index=False)
        print(filtered_df)
    else:
        merge_results.to_csv(Path(results_dir, 'statistics.csv'), index=False)
        print(merge_results)


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="statistical analysis of mobilome_tn_mining results")
    parser.add_argument('-d', "--results_dir", type=str, required=True, help="Path to the mobilome_tn_mining results directory")
    parser.add_argument("--min_copy", type=int, required=False, help="Minimum copy number filter")
    parser.add_argument("--max_copy", type=int, required=False, help="Maximum copy number filter")
    parser.add_argument("--min_orf", type=int, required=False, default=0, help="Minimal length of the ORF (nt)")
    parser.add_argument("--min_perc", type=float, required=False, default=0, help="Minimum percentage of the sequence that contains an ORF")

    # Parse arguments
    args = parser.parse_args()
    results_dir = Path(args.results_dir).resolve()

    # 统计log文件
    logs_folder = Path(results_dir, 'logs')
    analyze_logs_files(logs_folder,output_file=Path(results_dir, 'logstats.txt'))

    fasta_folder = Path(results_dir, 'fasta')

    if args.min_orf:
        orf_stats_df = run_ORFfinder(fasta_folder, results_dir, args)
        #合并统计结果
        #merge_results = pd.merge(fasta_df, orf_stats_df, on='filename')
        #merge_results['percent'] = merge_results[orfml] / merge_results['num_sequences']
        # 打印并保存结果
        #merge_results.to_csv(Path(results_dir, 'statistics.csv'), index=False)
        #print(merge_results)
    else:
        # 统计fasta文件
        fasta_df = analyze_fasta_files(fasta_folder,args.min_copy,args.max_copy)
        fasta_df.to_csv(Path(results_dir, 'statistics.csv'), index=False)
        print(fasta_df)

    # # 创建保存符合筛选条件序列的文件夹
    # fasta_filter_outdir = Path(results_dir, 'fasta_filter_out')
    # if fasta_filter_outdir.is_dir():
    #     shutil.rmtree(fasta_filter_outdir)
    # fasta_filter_outdir.mkdir(exist_ok=False, parents=False)

    # # 运行ORFfinder
    # orfml = "orf_gr_" + str(args.min_orf) + "nt"
    # orf_stats = []
    # for file in results_df['filename']:
    #     infile = Path(fasta_folder,file)
    #     outfile = Path(ORFfinder_outdir,file.replace('.fasta', '.orf.fasta'))
    #     cmd = f"ORFfinder -in {infile} -out {outfile} -ml {args.min_orf} -outfmt 0"
    #     run_cmd(cmd)
    #     if Path(outfile).stat().st_size == 0:
    #         outfile.unlink()
    #         continue
    #     filter_seq = Path(fasta_filter_outdir,file)
    #     sequences = []
    #     orf_records = list(SeqIO.parse(outfile, "fasta"))
    #     seqnames = {record.id.split("_")[1].split("(")[0] for record in orf_records}
    #     for record in SeqIO.parse(infile, "fasta"):
    #         if any(seqname in record.id for seqname in seqnames):
    #             sequences.append(record)
    #     SeqIO.write(sequences,filter_seq,"fasta")

    #     orf_stats.append({"filename":file, orfml:len(sequences)})



if __name__ == "__main__":
    main()