#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   mobilome_usearch_blastn_trim.py
@Time    :   2025/07/19 20:22:34
@Author  :   Naisu Yang 
@Version :   1.0
@Contact :   3298990@qq.com
'''

# here put the import lib

import os
import argparse
import subprocess
import shutil
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_uc_file(uc_file):
    """Parse UC file and return singleton and cluster records."""
    singletons_records = []
    clusters_records = []
    
    with open(uc_file, 'r') as f:
        for line in f:
            if line.startswith('C'):
                cluster_record, cluster_size = line.split('\t', 3)[1:3]
                if cluster_size == '1':
                    singletons_records.append(cluster_record)
                else:
                    clusters_records.append(cluster_record)
    return singletons_records, clusters_records

def merge_files(files, output_file):
    """Merge multiple files into one."""
    with open(output_file, 'wb') as out_f:
        for file in files:
            with open(file, 'rb') as in_f:
                shutil.copyfileobj(in_f, out_f)

def run_cmd(cmd):
    """Run shell command with error handling."""
    try:
        subprocess.run(cmd,stdout=None, stderr=None, check=True, shell=True)
    except Exception as e:
        print(f"Error running: {e}")
        raise

def run_mafft_parallel(fasta_dir, processes):

    mafft_out_dir = Path(fasta_dir).parent / f'mafft_{Path(fasta_dir).name}'
    mafft_out_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"find {fasta_dir} -maxdepth 1 -type f -name '*.fasta' | "
        f"parallel -j {processes} --bar "
        f"'mafft --auto --quiet {{}} > {mafft_out_dir}/$(basename {{}} .fasta).aln'"
    )

    print(cmd)
    print(f"run mafft...")

    run_cmd(cmd)

def process_clusters(clusters_records, cons, clusters, processes, flank_length, output):

    # Merge results
    cluster_blast_trim_merged_file = Path(output,'cluster_blast_trim.fasta')
    cluster_blast_trim_slop_merged_file = Path(output, f'cluster_blast_trim_slop{flank_length}.fasta')
    merge_files(Path(output,'clusters_trim_fasta').glob('*.fasta'), cluster_blast_trim_merged_file)
    print(f"Merged results saved to : {cluster_blast_trim_merged_file}")
    merge_files(Path(output,f'clusters_trim_slop{flank_length}_fasta').glob('*.fasta'), cluster_blast_trim_slop_merged_file)

    all_trim_merged_file = Path(output, f'all_trim.fasta')
    merge_files([cluster_blast_trim_slop_merged_file, singletons_trim_fasta_file], all_trim_merged_file)

def read_and_filter_consensus(cons, clusters_records, output_dir):
    """Read consensus sequences and filter based on cluster records."""
    cons_dict = {rec.id: rec for rec in SeqIO.parse(cons, 'fasta')}
    print(f"Total {len(cons_dict)} consensus seqs in {cons}")

    clusters_prefix_records = ['Cluster' + record for record in clusters_records]
    clusters_cons = {cons: record for cons, record in cons_dict.items() if cons in clusters_prefix_records}
    print(f"Found {len(clusters_cons)} clusters with >2 sequences")

    filtered_consensus_fasta = Path(output_dir, 'filtered_consensus.fasta')
    SeqIO.write(clusters_cons.values(), filtered_consensus_fasta, "fasta")
    print(f"Filtered consensus saved to {filtered_consensus_fasta}")
    
    return filtered_consensus_fasta

def split_consensus_files(clusters_cons, output_dir):
    """Split consensus files into individual files per cluster."""
    cluster_cons_split_dir = Path(output_dir, 'tmp', 'cluster_cons_split')
    cluster_cons_split_dir.mkdir(exist_ok=True)
    
    for rec in SeqIO.parse(clusters_cons, 'fasta'):
        output_path = os.path.join(cluster_cons_split_dir, rec.id[7:])
        with open(output_path, "w") as output_file:
            SeqIO.write(rec, output_file, "fasta")

    return cluster_cons_split_dir

def makeblastdb(fasta, blastdb_dir):
    """Create BLAST database from FASTA file."""
    db = Path(blastdb_dir,fasta.stem)
    cmd = (
        f"makeblastdb -in {fasta} -input_type fasta "
        f"-title {fasta.stem} -dbtype nucl -out {db} > /dev/null"
    )
    print(cmd)
    run_cmd(cmd)

    return db

def run_singleton_blast(filtered_consensus_fasta, singleton_file, output_dir, processes):

    blast_db_dir = Path(output_dir, 'tmp', 'blast_db')
    blast_db_dir.mkdir(parents=True, exist_ok=True) 
    cmd = (
        f"makeblastdb -in {singleton_file} -input_type fasta "
        f"-title $(basename {singleton_file} .fasta) -dbtype nucl -out {blast_db_dir}/$(basename {singleton_file} .fasta) -logfile /dev/null"
    )
    run_cmd(cmd)

    namemap_dir = Path(output_dir, 'tmp', 'namemap')
    namemap_dir.mkdir(parents=True, exist_ok=True)
    cmd = (
        f"grep \"^>\" {singleton_file} | awk -F\"[>() ]\" '{{print $2 \" \" $NF}}' > {namemap_dir}/singleton.namemap"
    )
    run_cmd(cmd)

    cmd = (
        f"blastn -query {filtered_consensus_fasta} -db {blast_db_dir}/$(basename {singleton_file} .fasta) "
        f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
        f"-task blastn -evalue 1e-5 -num_threads {processes} | "
        f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
        f"bedtools getfasta -fi {singleton_file} -bed /dev/stdin -fo /dev/stdout 2> /dev/null | "
        f"awk 'NR==FNR {{map[$1]=$2; next}} /^>/ {{key=substr($0,2); split(key, arr, \"(\"); "
        f"if (arr[1] in map) print $0 \" \" map[arr[1]]; else print $0; next}} {{print}}' "
        f"{namemap_dir}/singleton.namemap /dev/stdin > {Path(output_dir, 'singletons_trim.fasta')}"
    )

    print(f"Running BLASTN on singleton sequences...")

    run_cmd(cmd)

def run_clusters_blast(cluster_cons_split_dir, clusters_dir, output_dir, processes):

    tbl_dir = Path(output_dir, 'tmp', 'tbl')
    tbl_dir.mkdir(parents=True, exist_ok=True) 

    cmd = (
        f"find {cluster_cons_split_dir} -maxdepth 1 -type f | "
        f"parallel -j {processes} --bar "
        f"'blastn -query {{}} -subject {clusters_dir}/$(basename {{}}) "
        f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
        f"-task blastn -evalue 1e-5 > {tbl_dir}/$(basename {{}})'"
    )

    print(f"Running BLASTN on filtered consensus sequences...")
    run_cmd(cmd)

    return tbl_dir

def tbl2bed(tbl_dir, output_dir, processes):
    """Convert table format to BED format."""
    bed_dir = Path(output_dir, 'tmp', 'bed')
    bed_dir.mkdir(parents=True, exist_ok=True) 

    cmd = (
        f"find {tbl_dir} -maxdepth 1 -type f | "
        f"parallel -j {processes} "
        f"'mobilome_blast_tbl2bed -t {{}} -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 > {bed_dir}/$(basename {{}})'"
    )

    run_cmd(cmd)

    return bed_dir

def get_seq_len(clusters_dir, output_dir, processes):

    seq_len_dir = Path(output_dir, 'tmp', 'seq_len')
    seq_len_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"find {clusters_dir} -maxdepth 1 -type f -not -name '*.fai' | "
        f"parallel -j {processes} "
        f"'seqkit fx2tab -l -n -i {{}} > {seq_len_dir}/$(basename {{}})'"
    )

    run_cmd(cmd)

    return seq_len_dir

def bed_slop(bed_dir, seq_len_dir, flank_length, output_dir, processes):

    bed_slop_dir = Path(output_dir, 'tmp', f'bed_slop_{flank_length}')
    bed_slop_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"find {bed_dir} -maxdepth 1 -type f | "
        f"parallel -j {processes} "
        f"'bedtools slop -b {flank_length} -g {seq_len_dir}/$(basename {{}}) -i {{}} > {bed_slop_dir}/$(basename {{}})'"
    )

    run_cmd(cmd)

    return bed_slop_dir

def create_namemap(clusters_dir, output_dir, processes):
    namemap_dir = Path(output_dir, 'tmp', 'namemap')
    namemap_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"find {clusters_dir} -maxdepth 1 -type f -not -name '*.fai' | "
        f"parallel -j {processes} "
        f"'grep \"^>\" {{}} | awk -F\"[>() ]\" \"{{print \\$2 \\\" \\\" \\$NF}}\" > {namemap_dir}/$(basename {{}})'"
    )

    run_cmd(cmd)

    return namemap_dir

def get_clusters_trim_fasta(bed_dir, clusters_dir, namemap_dir, output_dir, processes):

    clusters_trim_fasta_dir = Path(output_dir, 'clusters_trim_fasta')
    clusters_trim_fasta_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"ls {bed_dir} | xargs -I {{}} echo "
        f"\"bedtools getfasta -bed {bed_dir}/{{}} -fi {clusters_dir}/{{}} | "
        f"awk 'NR==FNR {{map[\$1]=\$2; next}} /^>/ {{key=substr(\$0,2); split(key, arr, \\\"(\\\"); "
        f"if (arr[1] in map) print \$0 \\\" \\\" map[arr[1]]; else print \$0; next}} {{print}}' "
        f"{namemap_dir}/{{}} /dev/stdin > {clusters_trim_fasta_dir}/{{}}.fasta\" | parallel -j {processes} --bar"
    )
    print(f"get sequences...")
    run_cmd(cmd)

    return clusters_trim_fasta_dir

def get_clusters_trim_slop_fasta(bed_slop_dir, clusters_dir, namemap_dir, output_dir, processes):

    clusters_trim_slop_fasta_dir = Path(output_dir, 'clusters_trim_slop_fasta')
    clusters_trim_slop_fasta_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"ls {bed_slop_dir} | xargs -I {{}} echo "
        f"\"bedtools getfasta -bed {bed_slop_dir}/{{}} -fi {clusters_dir}/{{}} | "
        f"awk 'NR==FNR {{map[\$1]=\$2; next}} /^>/ {{key=substr(\$0,2); split(key, arr, \\\"(\\\"); "
        f"if (arr[1] in map) print \$0 \\\" \\\" map[arr[1]]; else print \$0; next}} {{print}}' "
        f"{namemap_dir}/{{}} /dev/stdin > {clusters_trim_slop_fasta_dir}/{{}}.fasta\" | parallel -j {processes} --bar"
    )
    print(f"get sequences...")
    run_cmd(cmd)

    return clusters_trim_slop_fasta_dir

def get_namemap(clusters_dir, output_dir):

    seq_ids_file = Path(output_dir, 'tmp', 'seq_ids.txt')
    seq_ids_namemap_file = Path(output_dir, 'tmp', 'seq_ids.namemap')

    seq_ids = {
        rec.id.partition(":")[0]
        for file in Path(clusters_dir).iterdir()
        if file.suffix != '.fai'
        for rec in SeqIO.parse(file, 'fasta')
    }

    seq_ids_file.write_text("\n".join(seq_ids))

    cmd = (
        f"rg -w -f {seq_ids_file} /storage/eukaryotes_genome/eukaryotes.namemap > {seq_ids_namemap_file}"
    )

    print('Search Genome accession...')

    run_cmd(cmd)

    with open('your_file.txt', 'r') as f:
        namemap = {line.strip().split()[1]: line.strip().split()[0] 
                for line in f if len(line.strip().split()) >= 2}

    return namemap

# def extend_clusters_flanking_regions(clusters_trim_fasta_dir, namemap, output_dir, processes):

#     extend_tmp_dir = Path(output_dir, 'tmp', 'extend_tmp')
#     extend_tmp_dir.mkdir(parents=True, exist_ok=True)

#     clusters_trim_extent_fasta_dir = Path(output_dir, 'clusters_trim_extent_fasta')
#     clusters_trim_extent_fasta_dir.mkdir(parents=True, exist_ok=True)

#     for file in Path(clusters_trim_fasta_dir).glob('*.fasta'):
#         for rec in SeqIO.parse(file, 'fasta'):

#     cmd = (
#         f"find {bed_dir} -maxdepth 1 -type f | "
#         f"parallel -j {processes} "
#         f"'bedtools slop -b 5000 -g  -i {{}} > {bed_slop_dir}/$(basename {{}})'"
#     )
    #-g /storage/eukaryotes_genome/fungi/GCA/GCA_025169535_1.fna.fai -i /dev/stdin | bedtools getfasta -fi /storage/eukaryotes_genome/fungi/GCA/GCA_025169535_1.fna -bed /dev/stdin -fo 3.fa

    # with open(file, 'rb') as in_f:
    #     shutil.copyfileobj(in_f, out_f)
    # for rec in SeqIO.parse(clusters_cons, 'fasta'):
    # output_path = os.path.join(cluster_cons_split_dir, rec.id[7:])
    # with open(output_path, "w") as output_file:
    #     SeqIO.write(rec, output_file, "fasta")
    # pass

def run_tirvish(cluster_cons_split_dir, clusters_dir, output_dir, processes):

    tirvish_tmp_dir = Path(output_dir, 'tmp', 'trivish_tmp')
    tirvish_tmp_dir.mkdir(parents=True, exist_ok=True)

    mafft_out_dir = Path(fasta_dir).parent / f'mafft_{Path(fasta_dir).name}'
    mafft_out_dir.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"find {fasta_dir} -maxdepth 1 -type f -name '*.fasta' | "
        f"parallel -j {processes} --bar "
        f"'mafft --auto --quiet {{}} > {mafft_out_dir}/$(basename {{}} .fasta).aln'"
    )

    print(cmd)
    print(f"run mafft...")

    run_cmd(cmd)

    mkdir -p tirvish/index && ls cluster-rename | xargs -I {} echo "gt suffixerator -db cluster-rename/{} -indexname tirvish/index/{} -tis -suf -lcp -des -ssp -sds -dna -mirrored" | parallel -j 40 --bar
    pass

def main():
    parser = argparse.ArgumentParser(
        description='Process USearch clustering results and perform BLASTN alignment.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    #parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-uc', required=True, help='USearch UC file')
    parser.add_argument('-cons', required=True, help='USearch Cluster consensus seqs')
    parser.add_argument('-clusters', required=True, help='USearch clusters directory')
    parser.add_argument('--flank_length', type=int, default=50, help="Length to extend on both sides of the DNA sequence (default: 50).")
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-p', '--processes', type=int, default=1, help='number of parallel processes to use')
    
    args = parser.parse_args()
    
    # Setup output directory
    output_dir = Path(args.output)
    if output_dir.exists() and output_dir.is_dir():
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)
    Path(output_dir, 'tmp').mkdir(exist_ok=True)

    # Parse UC file
    singletons_records, clusters_records = parse_uc_file(args.uc)
    print(f"Found {len(singletons_records)} singletons and {len(clusters_records)} clusters in {args.uc}")

    #Process singletons and clusters
    singleton_file = Path(output_dir, 'tmp', 'singletons.fasta')
    files = [Path(args.clusters,record) for record in singletons_records]
    merge_files(files, singleton_file)
    print(f"Saved {len(singletons_records)} singletons to {singleton_file}")

    clusters_file = Path(output_dir, 'tmp', 'clusters.fasta')
    files = [Path(args.clusters,record) for record in clusters_records]
    merge_files(files, clusters_file)
    print(f"Saved {len(clusters_records)} clusters to {clusters_file}")

    get_namemap(args.clusters, output_dir,  args.processes)

    filtered_consensus_fasta = read_and_filter_consensus(args.cons, clusters_records, output_dir)

    cluster_cons_split_dir = split_consensus_files(filtered_consensus_fasta, output_dir)
    
    run_singleton_blast(filtered_consensus_fasta, singleton_file, output_dir, args.processes)

    tbl_dir = run_clusters_blast(cluster_cons_split_dir, args.clusters, output_dir, args.processes)

    bed_dir = tbl2bed(tbl_dir, output_dir, args.processes)

    namemap_dir = create_namemap(args.clusters, output_dir, args.processes)

    seq_len_dir = get_seq_len(args.clusters, output_dir, args.processes)

    bed_slop_dir = bed_slop(bed_dir, seq_len_dir, args.flank_length, output_dir, args.processes)

    clusters_trim_fasta_dir = get_clusters_trim_fasta(bed_dir, args.clusters, namemap_dir, output_dir, args.processes)


    # namemap = get_namemap(args.clusters, output_dir)
    # extend_clusters_flanking_regions(clusters_trim_fasta_dir, namemap, output_dir, args.processes)
    #run_mafft_parallel(clusters_trim_fasta_dir, args.processes)
    #clusters_trim_slop_fasta_dir = get_clusters_trim_slop_fasta(bed_slop_dir, args.clusters, namemap_dir, output_dir, args.processes)
    #run_mafft_parallel(clusters_trim_slop_fasta_dir, args.processes)


if __name__ == '__main__':
    main()