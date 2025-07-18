#!/usr/bin/env python3
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

def tbl2bed(tbl,bed):
    """Convert table format to BED format."""
    with open(tbl) as f1,open(bed,'w') as f2:
        for line in f1:
            if not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 10:
                    chrom = parts[1]
                    chrom_start = int(parts[8])
                    chrom_end = int(parts[9])
                if chrom_start > chrom_end:
                    line = f'{chrom}\t{chrom_end - 1}\t{chrom_start}\t.\t.\t-\n'
                else:
                    line = f'{chrom}\t{chrom_start - 1}\t{chrom_end}\t.\t.\t+\n'
                f2.write(line)

def run_cmd(cmd):
    """Run shell command with error handling."""
    try:
        subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True,shell=True)
    except Exception as e:
        print(f"Error running: {e}")
        raise


def run_blastn(query, subject, output, flank_length):

    """Run BLASTN with pipeline processing."""

    if flank_length:
        bed_dir = Path(output) / 'tmp' / 'bed'
        bed_dir.mkdir(exist_ok=True)
        bed_file = (bed_dir / query.stem).with_suffix('.bed')
        bed_slop_file = (bed_dir / query.stem).with_suffix(f'.slop{flank_length}.bed')
        len_file = (bed_dir / subject.stem).with_suffix('.len')
        cmd = (
        f"seqkit fx2tab -l -n -i {subject} > {len_file}"
        )
        run_cmd(cmd)

        cmd = (
        f"blastn -query {query} -subject {subject} -strand plus "
        f"-outfmt 6 -max_target_seqs 99999999 -task blastn -evalue 1e-5 | "
        f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 > {bed_file} "
    )
        run_cmd(cmd)

        cmd = (
        f"bedtools slop -s -l {flank_length} -r {flank_length} -g {len_file} -i {bed_file} > {bed_slop_file} "
   
    )
        run_cmd(cmd)

        blast_clusters_out_dir = Path(output,'clusters_trim_fasta')
        blast_clusters_slop_out_dir = Path(output, f'clusters_trim_slop{flank_length}_fasta')
        blast_clusters_output_file = (blast_clusters_out_dir / query.stem).with_suffix('.fasta')
        blast_clusters_slop_output_file = (blast_clusters_slop_out_dir / query.stem).with_suffix(f'.slop{flank_length}.fasta')

        cmd = (
            f"bedtools getfasta -fi {subject} -bed {bed_file} -fo {blast_clusters_output_file}"
        )
        run_cmd(cmd)

        cmd = (
            f"bedtools getfasta -fi {subject} -bed {bed_slop_file} -fo {blast_clusters_slop_output_file}"
        )
        run_cmd(cmd)

    # cmd = (
    #     f"blastn -query {query} -subject {subject} -strand plus "
    #     f"-outfmt 6 -max_target_seqs 99999999 -task blastn -evalue 1e-5 | "
    #     f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
    #     f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
    #     f"bedtools getfasta -fi {subject} -bed /dev/stdin -fo {output_file}"
    # )
    # run_cmd(cmd)

def makeblastdb(fasta, blastdb_dir):
    """Create BLAST database from FASTA file."""
    db = Path(blastdb_dir,fasta.stem)
    cmd = (
        f"makeblastdb -in {fasta} -input_type fasta "
        f"-title {fasta.stem} -dbtype nucl -out {db}"
    )
    run_cmd(cmd)

    return db

def run_mafft(fasta, out):
    cmd = f"mafft --auto {fasta} > {out}"
    run_cmd(cmd)

def run_singleton_blast(query, singleton_file, db, output_file, processes):

    cmd = (
        f"blastn -query {query} -db {db} "
        f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
        f"-task blastn -evalue 1e-5 -num_threads 1 | "
        f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
        f"bedtools getfasta -fi {singleton_file} -bed /dev/stdin -fo {output_file}"
    )
    run_cmd(cmd)

def process_clusters(clusters_records, cons, clusters, processes, flank_length, output):

    # Setup directories
    tmp_dir = Path(output,'tmp')
    tmp_dir.mkdir(exist_ok=True)

    # Read consensus sequences
    cons_dict = {rec.id: rec for rec in SeqIO.parse(cons, 'fasta')}
    print(f"Total {len(cons_dict)} consensus seqs in {cons}")

    # Filter cluster consensus sequences
    clusters_prefix_records = ['Cluster' + record for record in clusters_records]
    clusters_cons = {cons: record for cons, record in cons_dict.items() if cons in clusters_prefix_records}
    print(f"Found {len(clusters_cons)} clusters with >2 sequences")

    # Write filtered consensus
    filtered_consensus_fasta = Path(output,'filtered_consensus.fasta')
    SeqIO.write(clusters_cons.values(), filtered_consensus_fasta, "fasta")
    print(f"Filtered consensus saved to {filtered_consensus_fasta}")

    # Split consensus files
    cluster_cons_split_dir = Path(tmp_dir,'cluster_cons_split')
    cluster_cons_split_dir.mkdir(exist_ok=True)
    for cluster, record in clusters_cons.items():
        output_path = os.path.join(cluster_cons_split_dir, cluster)
        with open(output_path, "w") as output_file:
            SeqIO.write(record, output_file, "fasta")

    #makeblastdb for singleton
    blastdb_dir = Path(tmp_dir,'db')
    blastdb_dir.mkdir(exist_ok=True)
    singleton_file = Path(tmp_dir, 'singletons.fasta')
    db = makeblastdb(singleton_file,blastdb_dir)

    # Prepare parallel processing
    blast_clusters_args_list = []
    blast_singletons_args_list = []
    blast_clusters_mafft_args_list = []
    blast_clusters_slop_mafft_args_list = []
    blast_singletons_mafft_args_list = []

    blast_clusters_out_dir = Path(output,'clusters_trim_fasta')
    blast_clusters_slop_out_dir = Path(output, f'clusters_trim_slop{flank_length}_fasta')
    #blast_singletons_out_dir = Path(output,'singletons_trim_fasta')
    blast_clusters_mafft_out_dir = Path(output,'clusters_trim_fasta_mafft_out')
    blast_clusters_slop_mafft_out_dir = Path(output,f'clusters_trim_slop{flank_length}_fasta_mafft_out')
    #blast_singletons_mafft_out_dir = Path(output,'singletons_trim_fasta_mafft_out')

    blast_clusters_out_dir.mkdir(exist_ok=True)
    blast_clusters_slop_out_dir.mkdir(exist_ok=True)
    #blast_singletons_out_dir.mkdir(exist_ok=True)
    blast_clusters_mafft_out_dir.mkdir(exist_ok=True)
    blast_clusters_slop_mafft_out_dir.mkdir(exist_ok=True)
    #blast_singletons_mafft_out_dir.mkdir(exist_ok=True)

    for cluster_record in clusters_records:
        query = Path(cluster_cons_split_dir, 'Cluster' + cluster_record)
        subject = Path(clusters, cluster_record)

        blast_clusters_output_file = Path(blast_clusters_out_dir , 'Cluster' + cluster_record + ".fasta")
        blast_clusters_slop_output_file = Path(blast_clusters_slop_out_dir , 'Cluster' + cluster_record + f".slop{flank_length}.fasta")
        #blast_singletons_output_file = Path(blast_singletons_out_dir , 'Cluster' + cluster_record + ".fasta")
        blast_clusters_mafft_output_file = Path(blast_clusters_mafft_out_dir , 'Cluster' + cluster_record + ".aln")
        blast_clusters_slop_mafft_output_file = Path(blast_clusters_slop_mafft_out_dir , 'Cluster' + cluster_record + f".slop{flank_length}.aln")
        
        blast_clusters_args_list.append((query, subject, output,flank_length))
        #blast_singletons_args_list.append((query, singleton_file, db, blast_singletons_output_file, processes))
        blast_clusters_mafft_args_list.append((blast_clusters_output_file, blast_clusters_mafft_output_file))
        blast_clusters_slop_mafft_args_list.append((blast_clusters_slop_output_file, blast_clusters_slop_mafft_output_file))


    # Run BLASTN in parallel
    print(f"Running BLASTN on filtered consensus sequences...")
    with Pool(processes=processes) as pool:
        pool.starmap(run_blastn, blast_clusters_args_list)
    print(f"BLASTN completed. Output in: {blast_clusters_out_dir}")

    # Merge results
    cluster_blast_trim_merged_file = Path(output,'cluster_blast_trim.fasta')
    cluster_blast_trim_slop_merged_file = Path(output, f'cluster_blast_trim_slop{flank_length}.fasta')
    merge_files(Path(output,'clusters_trim_fasta').glob('*.fasta'), cluster_blast_trim_merged_file)
    print(f"Merged results saved to : {cluster_blast_trim_merged_file}")
    merge_files(Path(output,f'clusters_trim_slop{flank_length}_fasta').glob('*.fasta'), cluster_blast_trim_slop_merged_file)

    # BLASTN on singletons
    print(f"Running BLASTN on singleton sequences...")

    # cmd = (
    #     f"ls {cluster_cons_split_dir} | xargs -I {{}} echo '"
    #     f"blastn -query {cluster_cons_split_dir}/{{}} -db {db} "
    #     f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
    #     f"-task blastn -evalue 1e-5 -num_threads 1 | "
    #     f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
    #     f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
    #     f"bedtools getfasta -fi {singleton_file} -bed /dev/stdin -fo {blast_singletons_out_dir}/{{}}' | parallel -j {processes} --bar  "
    # )
    singletons_trim_fasta_file = Path(output,'singletons_trim_fasta.fasta')
    cmd = (
        f"blastn -query {filtered_consensus_fasta} -db {db} "
        f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
        f"-task blastn -evalue 1e-20 -num_threads {processes} | "
        f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
        f"bedtools getfasta -fi {singleton_file} -bed /dev/stdin -fo {singletons_trim_fasta_file}"
    )
    try:
        subprocess.run(cmd,stdout=None,stderr=None,check=True,shell=True)
    except Exception as e:
        print(f"Error running: {e}")
        raise

    all_trim_merged_file = Path(output, f'all_trim.fasta')
    merge_files([cluster_blast_trim_slop_merged_file, singletons_trim_fasta_file], all_trim_merged_file)

    # with Pool(processes=processes) as pool:
    #     pool.starmap(run_singleton_blast, blast_singletons_args_list)
    # print(f"BLASTN completed. Output in: {blast_singletons_out_dir}")

    # valid_files = []
    # # 遍历文件夹中的所有文件
    # for file in Path(blast_singletons_out_dir).iterdir():
    #     if file.is_file():
    #         if file.stat().st_size == 0:
    #             # 删除大小为0的文件
    #             file.unlink()
    #         else:
    #             valid_files.append(file)
    # for file in valid_files:
    #     blast_singletons_mafft_output_file = Path(blast_singletons_mafft_out_dir , file.stem + ".aln")
    #     blast_singletons_mafft_args_list.append((file, blast_singletons_mafft_output_file))

    # Run MAFFT in parallel
    print(f"Running MAFFT on each cluster...")
    with Pool(processes=processes) as pool:
        pool.starmap(run_mafft, blast_clusters_mafft_args_list)
        pool.starmap(run_mafft, blast_clusters_slop_mafft_args_list)
        #pool.starmap(run_mafft, blast_singletons_mafft_args_list)
    print(f"MAFFT completed. Output in: {blast_clusters_mafft_out_dir}")

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
        output_path = os.path.join(cluster_cons_split_dir, rec.id)
        with open(output_path, "w") as output_file:
            SeqIO.write(rec, output_file, "fasta")

    return cluster_cons_split_dir

def run_clusters_blast(cluster_cons_split_dir, clusters_dir, output_file, processes):

    cmd = (
        f"blastn -query {query} -db {db} "
        f"-strand plus -outfmt 6 -max_target_seqs 99999999 "
        f"-task blastn -evalue 1e-5 -num_threads 1 | "
        f"mobilome_blast_tbl2bed -t /dev/stdin -o /dev/stdout | "
        f"sort -k1,1 -k2,2n | bedtools merge -d 10000 | "
        f"bedtools getfasta -fi {singleton_file} -bed /dev/stdin -fo {output_file}"
    )
    run_cmd(cmd)

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

    filtered_consensus_fasta = read_and_filter_consensus(args.cons, clusters_records, output_dir)

    split_consensus_files(filtered_consensus_fasta, output_dir)
    # Process clusters with BLASTN and MAFFT
    #process_clusters(clusters_records, args.cons, args.clusters, args.processes, args.flank_length, args.output)

if __name__ == '__main__':
    main()