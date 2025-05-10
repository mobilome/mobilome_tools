#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   mobilome_create_tn_mining_jobs.py
@Time    :   2025/05/03 17:28:37
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com
'''

# here put the import lib

import os
import argparse
from pathlib import Path

def count_lines(filename):
    with open(filename, 'r') as file:
        count = 0
        for _ in file:
            if _.strip():
                count += 1
    return count

def create_slurm_script(args,path):
    GCA_list = Path(path,'GCA.list')
    job_number = count_lines(GCA_list)
    slurm_script_content = f"""#!/bin/bash
#SBATCH --job-name={path.name}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/dev/null
#SBATCH --array=1-{job_number}%48

PERFIX="{args.prefix}"
CATEGORY="{path.parent.name}/{path.name}"

#输出目录，默认家目录下的PERFIX
OUTPUT="$HOME/$PERFIX/$CATEGORY"
FASTADIR="$OUTPUT/fasta"
LOGDIR="$OUTPUT/logs"
TMPDIR="$OUTPUT/tmpdir"
#创建目录
mkdir -p "$OUTPUT" "$FASTADIR" "$LOGDIR" "$TMPDIR"

# 获取任务ID
TASK_ID=$SLURM_ARRAY_TASK_ID

# 获取基因组文件
mapfile -t GENOME_LIST < {GCA_list}
GENOME_FILE=${{GENOME_LIST[$TASK_ID - 1]}}

#检查文件是否存在
if [ ! -f "$GENOME_FILE" ]; then
   echo "Error: Genome file $GENOME_FILE not found"
   exit 1
fi

#运行mobilome_tn_mining
echo mobilome_tn_mining -query {args.dde} -db ${{GENOME_FILE%.fna}} -name_id {Path(path,'name2id.txt')}  -tmp_dir $TMPDIR -left_len {args.left_len} -right_len {args.right_len} -out $FASTADIR/$(basename "$GENOME_FILE" .fna).fasta -num 8 -tblastn_option {args.blast_options} > $LOGDIR/$(basename "$GENOME_FILE" .fna).log

mobilome_tn_mining -query {args.dde} -db ${{GENOME_FILE%.fna}} -name_id {Path(path,'name2id.txt')}  -tmp_dir $TMPDIR -left_len {args.left_len} -right_len {args.right_len} -out $FASTADIR/$(basename "$GENOME_FILE" .fna).fasta -num 8 -tblastn_option {args.blast_options} >> $LOGDIR/$(basename "$GENOME_FILE" .fna).log
"""

    # 将内容写入文件
    slurm_job_file = Path(Path(args.outdir),f'{path.name}.slurm')
    with open(slurm_job_file, 'w') as f:
        f.write(slurm_script_content)

    print(f"总任务数: {job_number}")
    print(f"SLURM任务文件位置: {slurm_job_file}")
    print(f"DDE文件: {args.dde}")
    print(f"左侧延伸长度: {args.left_len},右侧延伸长度: {args.right_len}")
    print(f"blastn参数: {args.blast_options}")
    print(f"结果输出目录: {Path.home()}/{args.prefix}/{path.parent.name}/{path.name}")
    print("请检查脚本内容，然后使用以下命令提交作业:")
    print(f"sbatch {slurm_job_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate SLURM submission script for mobilome_tn_mining pipeline')

    parser.add_argument(
        '--prefix',
        required=True,
        help='Project prefix, suggest using tn name (e.g., "academ")'
    )
    parser.add_argument(
        '--dde',
        required=True,
        help='Path to DDE protein sequence file (FASTA format)'
    )
    parser.add_argument(
        '--left_len',
        type=int,
        default=3000,
        help='Number of base pairs to extend upstream (left)'
    )
    parser.add_argument(
        '--right_len',
        type=int,
        default=3000,
        help='Number of base pairs to extend downstream (right)'
    )
    parser.add_argument(
        '--blast_options',
        default='-evalue 1e-20',
        help='BLAST search options'
    )
    parser.add_argument(
        '--outdir',
        default='.',
        help='Output file directory'
    )

    args = parser.parse_args()

    Path(args.outdir).mkdir(exist_ok=True, parents=True)
    directories = sorted([child for child in Path("/home/bpool/storage/eukaryotes_genome").iterdir() if child.is_dir()])

    print(f"\n{'*' * 20}共生成 {len(directories)} 个SLURM任务文件{'*' * 20}")


    for index, directory in enumerate(directories, start=1):
        print(f"\n{'=' * 20}第{index}个SLURM任务文件信息{'=' * 20}")
        create_slurm_script(args, Path(directory))


if __name__ == "__main__":
    main()

