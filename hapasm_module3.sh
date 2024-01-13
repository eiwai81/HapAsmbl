#!/usr/bin/env bash

# Usage:
#    parallel -j 2 " bash ./scripts/hapasm_module3.sh {1} {2} {3} {4} \
#       haplotype_assembly/{1}/{2}/assm \
#       haplotype_assembly/{1}/{2}/{2}.{1}.{3}.fastq.gz" ::: VRN1 ::: barcode01 ::: h1 h2 ::: 2

# Dead programs tell no lies
set -ueo pipefail

# Get commandline arguments
region=$1
barcode=$2
hap_num=$3
algorithm=$4
asm_dir=$5
input_fq=$6
# set environment variables

# create directory for assembly
mkdir -p ${asm_dir} 

# Prefix of output consensus sequences
prefix=$(echo "${barcode}.${region}.${hap_num}")

# 1. USE SPOA TO GENERATE CONSENSUS SEQUENCE(S)
# spoa --strand-ambiguous --algorithm 1 reads.fastq > out.fasta
# -s,  --strand-ambiguous: for each sequence pick the strand with the better alignment
# -l, --algorithm:  2 - semi-global

echo "Generate consensus sequence from ${input_fq}"
spoa -s -l ${algorithm} ${input_fq} > ${asm_dir}/${prefix}.con.fasta

# 2. POLISH SPOA CONSENSUS WITH FLYE POLISHER
# flye --polish-target SEQ_TO_POLISH --pacbio-raw READS --iterations NUM_ITER --out-dir OUTPUTDIR --threads THREADS
# I used 5 iterations (--iterations 5)

echo "Polishing ${asm_dir}/${prefix}.con.fasta"
flye --polish-target ${asm_dir}/${prefix}.con.fasta --nano-raw ${input_fq} --iterations 5 --out-dir ${asm_dir}/${hap_num}_flye_polish

# 3. PROCESSING AFTER ASSEMBLY AND POLISHING
# a.Rename sequence header in flye-polished sequences from consensus to barcode.region.hap#.p.fasta

echo "Adding sample haplotype information to polished sequences"
cat ${asm_dir}/${hap_num}_flye_polish/polished_5.fasta | sed -r -e "s@Consensus@${prefix}@g" > ${asm_dir}/${prefix}.fasta

rm -rf ${asm_dir}/${hap_num}_flye_polish
rm -rf ${asm_dir}/${prefix}.con.fasta*

echo "Done."
echo "Done.."