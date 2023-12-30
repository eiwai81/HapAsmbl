#!/usr/bin/env python3
"""
Author : Ayodele Fakoya <eiwai@localhost>
Date   : 2023-06-14
Purpose: Fetch reads of homozygous and heterozygous samples
"""
import os
import re
import subprocess
import pandas as pd
import argparse
from typing import NamedTuple
from functools import partial


class Args(NamedTuple):
    """ Command-line arguments """
    input_tsv: str
    region: str
    barcode_id: str
    haplotagged_bam_dir: str
    fq_dir: str

# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Fetch reads of homozygous and heterozygous samples',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_tsv',
                        metavar='haplotag_tsv_file',
                        help='tsv containing haplotag information of reads')

    parser.add_argument('-r',
                        '--region',
                        help='target region',
                        metavar='region',
                        type=str,
                        default='none')

    parser.add_argument('-b',
                        '--barcode_id',
                        help='barcode identifier',
                        metavar='barcode',
                        type=str,
                        default='none')

    parser.add_argument('-o',
                        '--haplotagged_bam_dir',
                        help='Path of haplotagged bamfile',
                        metavar='BAM DIR',
                        type=str,
                        default='none')

    parser.add_argument('-f',
                        '--fq_dir',
                        help='Path to sample FASTQ file',
                        metavar='FASTQ DIR',
                        type=str,
                        default='none')

    args = parser.parse_args()

    return Args(args.input_tsv,
                args.region,
                args.barcode_id,
                args.haplotagged_bam_dir,
                args.fq_dir)


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()

    print(f'input_tsv: {args.input_tsv}')
    print(f'region: {args.region}')
    print(f'barcode_id: {args.barcode_id}')
    
    col_list, elem_count = parse_tsv(args.input_tsv)
    
    # Case 1: HOM_REF
    if elem_count == 0:
        print(f'{args.input_tsv} is HOM_REF')

        # Path to haplotype assembly and haplotagged bams
        bc_bam_dir = os.path.join(args.haplotagged_bam_dir,
                                  args.region,
                                  args.barcode_id)
        in_bam_list = list(
            filter(partial(re.match, r'barcode.*haplotagged\.bam$'),
                   os.listdir(f'{bc_bam_dir}')))

        # 1. Extract read ids
        for file in in_bam_list:
            assert os.path.isfile(f'{bc_bam_dir}/{file}'), 'FILE NOT FOUND'
            print('bam_file:', os.path.basename(f'{bc_bam_dir}/{file}'))

            samtools_cmd=f"""
            samtools view {bc_bam_dir}/{file} \
            | cut -f 1 \
            > {bc_bam_dir}/read_ids.txt \
            """.strip('\n').replace('    ','')
            subprocess.run(samtools_cmd, shell=True, check=True)

        # Extract read with seqkit
        ## Path to sample fastq files
        bc_fq_dir = os.path.join(os.getcwd(), args.fq_dir)
        in_fq_list=list(
            filter(partial(re.match, r'barcode.*filt\.fastq'), 
                   os.listdir(args.fq_dir)))

        for file in in_fq_list:
            if file.startswith(args.barcode_id):
                assert os.path.isfile(f'{bc_fq_dir}/{file}')
                print('fastq_file:', os.path.basename(f'{bc_fq_dir}/{file}'))

                get_h1_fqs=f"""
                seqkit grep --pattern-file {bc_bam_dir}/read_ids.txt \
                {bc_fq_dir}/{file} \
                > {bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz
                """.strip().replace('    ', '')
                subprocess.run(get_h1_fqs, shell=True, check=True)
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz')
            
                get_h2_fqs=f"""
                seqkit grep --pattern-file {bc_bam_dir}/read_ids.txt \
                {bc_fq_dir}/{file} \
                > {bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz
                """.strip().replace('    ', '')
                subprocess.run(get_h2_fqs, shell=True, check=True)
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz')
        
    # Case 2 - HOM_ALT
    elif elem_count != 0 and ("H1" not in col_list) and ("none" in col_list):
        print(f'{args.input_tsv} is HOM_ALT')

        # Path to haplotype assembly and haplotagged bams
        bc_bam_dir = os.path.join(args.haplotagged_bam_dir,
                                  args.region,
                                  args.barcode_id)
        in_bam_list = list(
            filter(partial(re.match, r'barcode.*haplotagged\.bam$'),
                   os.listdir(f'{bc_bam_dir}')))

        # 1. Extract read ids
        for file in in_bam_list:
            assert os.path.isfile(f'{bc_bam_dir}/{file}'), 'FILE NOT FOUND'
            print('bam_file:', os.path.basename(f'{bc_bam_dir}/{file}'))

            samtools_cmd=f"""
            samtools view {bc_bam_dir}/{file} \
            | cut -f 1 \
            > {bc_bam_dir}/read_ids.txt \
            """.strip('\n').replace('    ','')
            subprocess.run(samtools_cmd, shell=True, check=True)

        # Extract read with seqkit
        ## Path to sample fastq files
        bc_fq_dir = os.path.join(os.getcwd(), args.fq_dir)
        in_fq_list=list(
            filter(partial(re.match, r'barcode.*filt\.fastq'), 
                   os.listdir(args.fq_dir)))

        for file in in_fq_list:
            if file.startswith(args.barcode_id):
                assert os.path.isfile(f'{bc_fq_dir}/{file}')
                print('fastq_file:', os.path.basename(f'{bc_fq_dir}/{file}'))

                get_h1_fqs=f"""
                seqkit grep --pattern-file {bc_bam_dir}/read_ids.txt \
                {bc_fq_dir}/{file} \
                > {bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz
                """.strip().replace('    ', '')
                subprocess.run(get_h1_fqs, shell=True, check=True)
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz')
            
                get_h2_fqs=f"""
                seqkit grep --pattern-file {bc_bam_dir}/read_ids.txt \
                {bc_fq_dir}/{file} \
                > {bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz
                """.strip().replace('    ', '')
                subprocess.run(get_h2_fqs, shell=True, check=True)
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz')
    
    # Case 3: HET
    elif elem_count != 0 and ("H1" in col_list) and ("H2" in col_list):
        print(f'{args.input_tsv} is HET')
        
        # Path to haplotype assembly and haplotagged bams
        bc_bam_dir = os.path.join(args.haplotagged_bam_dir,
                                  args.region,
                                  args.barcode_id)
        bc_fq_dir = os.path.join(os.getcwd(), args.fq_dir)
        in_fq_list=list(
            filter(partial(re.match, r'barcode.*filt\.fastq'), 
                   os.listdir(args.fq_dir)))
        
        for file in in_fq_list:
            if file.startswith(args.barcode_id):
                assert os.path.isfile(f'{bc_fq_dir}/{file}')
                print('fastq_file:', os.path.basename(f'{bc_fq_dir}/{file}'))

                whatshap_split_cmd=f"""
                whatshap split \
                --output-h1 {bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz \
                --output-h2 {bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz \
                {bc_fq_dir}/{file} \
                {args.input_tsv} 
                """.strip('\n').replace('    ','')
                subprocess.run(whatshap_split_cmd, shell=True, check=True)
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h1.fastq.gz')
                assert os.path.isfile(f'{bc_bam_dir}/{args.barcode_id}.{args.region}.h2.fastq.gz')

# --------------------------------------------------
def parse_tsv(input_tsv: str):
    """ Parse TSV file"""
    df=pd.read_csv(input_tsv, sep='\t',
                dtype={"haplotype":"string",
                       "phaseset":"string",
                       "chromosome":"string"})
    
    col_list=list(df["haplotype"])
    
    elem_count=len(col_list)
    
    return col_list, elem_count

# --------------------------------------------------
if __name__ == '__main__':
    main()
