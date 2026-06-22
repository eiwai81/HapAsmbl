#!/usr/bin/env bash

# dead programs tell no lies
set -ueo pipefail

# Get options
# Usage information
# check whether user had supplied -h or --help . If yes display usage
if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
	echo "Usage: $0 [-r reference_file] [-f fastq_file] [-o output_dir] [-b barcode_id] [-t threads] [-m flowcell_chemistry]"
	exit 0
fi

# Get commandline arguments
while getopts r:f:o:b:t:m: flag
do
    case "${flag}" in
        r) reference=${OPTARG};;
        f) fastq_file=${OPTARG};;
        o) out_dir=${OPTARG};;
        b) barcode_id=${OPTARG};;
        t) threads=${OPTARG};;
        m) model_config=${OPTARG};;
    esac
done

# Set current working directory
working_dir=$(echo "$(pwd)")

# Set output directory
hapasm_dir=$(echo "${out_dir}")

# print genes (or sequence IDs) in the reference and save to a text file
seqkit seq --name --only-id ${reference} > ${working_dir}/seq_ids.txt

# path to file containing names of sequences in reference file
seq_ids=$(echo "${working_dir}/seq_ids.txt")

# Check if $seq_ids exists and is not an empty file,
if [ -f "$seq_ids" ] && [ -s "$seq_ids" ]; then
    echo "seq_ids.txt is not empty."
else
    echo "Error: seq_ids.txt cannot be empty."
    exit 1
fi

# Create directories for sample bam and vcf files
mkdir -p ${hapasm_dir}/{bam,vcf}
mkdir -p ${hapasm_dir}/bam/${barcode_id}

######## MODULE 1 ##########################

# ---------- 1.2 Map reads to reference ---------------
echo "Mapping filtered (unmapped) reads to the reference file"
minimap2 --MD -a -x map-ont ${reference} ${fastq_file} | samtools sort > ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_tmp.bam

# Sort indexed bam file
samtools index ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_tmp.bam

echo "Removing unmapped and chimeric reads"
# Filter out umapped reads with -F 2308 flag 
# Include -h flag in samtools view to include header information 
## Exclude chimeric reads with samtools with bbmap's reformat.sh
## clipfilter=10 discards reads with more than 10 soft-clipped bases

samtools view -h -F 2308 ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_tmp.bam \
| reformat.sh clipfilter=10 in=stdin.bam out=stdout.bam \
| samtools sort > ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam

samtools index ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam

# remove unfiltered bam file
rm -rf ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_tmp.bam*

######## MODULE 2 ##########################

# ---------- 2.1 Call variants - clair3  ---------------
# Call variants (SNPs and small indels) with clair3
platform="ont"

# Set model for variant calling
if [[ "${model_config}" == "R9.4.1" ]]; then
	# Set model for R9.4.1 data
    model_path=$(echo "${CONDA_PREFIX}/bin/models/r941_prom_sup_g5014")
elif [[ "${model_config}" == "R10.4.1" ]]; then
	# Set model for R10.4.1 data
    model_path=$(echo "$CONDA_PREFIX/bin/models/r1041_e82_400bps_sup_v500")
else
    echo "Error: Invalid flow cell chemistry. Please use -m R9.4.1 for R9.4.1 LRAS data or -m R10.4.1 for R10.4.1 LRAS data."
    exit 1
fi

# path to bed
bed_path=$(echo "$(pwd)/ref")

# Modification
## I removed '--bed_fn=fn.bed' and replaced with --ctg_name=STR to
## call variants along the whole gene
run_clair3.sh \
--bam_fn=${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam \
--ref_fn=${reference} \
--threads=${threads} \
--platform=${platform} \
--model_path=${model_path} \
--output=${hapasm_dir}/vcf/${barcode_id} \
--include_all_ctgs \
--sample_name=${barcode_id} \
--chunk_size=25000 \
--var_pct_full=1 \
--ref_pct_full=1 \
--snp_min_af=0.01 \
--no_phasing_for_fa \
--use_whatshap_for_final_output_phasing \
--enable_variant_calling_at_sequence_head_and_tail \
--remove_intermediate_dir

# ---------- 3.1 Phase variants - whatshap ---------------
whatshap phase \
-o ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_tmp.vcf.gz \
--reference ${reference} \
--tag HP \
${hapasm_dir}/vcf/${barcode_id}/merge_output.vcf.gz \
${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam \
--indels \
--sample ${barcode_id} \
--ignore-read-groups \
--internal-downsampling 23 \
--distrust-genotypes

# FIX: Remove variants in heterozygous where HP tage is missing, which leads to NoneType Error in haplotag step
bcftools view -i 'GT="het" && FMT/HP!="."' ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_tmp.vcf.gz \
-Oz -o ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_phased.vcf.gz

# Index phased compressed vcf.gz file
tabix -f -p vcf ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_phased.vcf.gz

# Housekeeping
rm -rf ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_tmp.vcf.gz

#------- 3.2 Haplotag reads in bamfile using phased VCF and ----------------------------
#-------     generate list of reads belonging to each haplotype group
whatshap haplotag \
-o ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotagged.bam \
--reference ${reference} \
--output-haplotag-list ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotags.tsv \
--ignore-read-groups \
--sample ${barcode_id} \
--skip-missing-contigs \
${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_phased.vcf.gz \
${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam

# Index haplotagged bamfile
samtools index ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotagged.bam

#-------- 3.3 Per gene, cluster reads from each haplotype ---------
# make directory for each gene with sub-directory for each sample
cat ${seq_ids} | parallel -j 1 "mkdir -p ${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp"

# Get haplotagged reads per gene
my_awk=$(echo 'BEGIN {OFS = "\t"} /^#/ {print} !/^#/ && $4 ~ var {print}')

cat ${seq_ids} | parallel "cat ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotags.tsv \
| awk -v var={1} '$my_awk' > ${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp/${barcode_id}_{1}.haplotags.tsv" 

# -------- 4 Per sample, split reads originating from each flowering gene in haplotagged bamfile -------
# 1. per sample, extract reads that map to each gene region listed in the seq_ids.txt file
# use -h flag in samtools view to include header
cat ${seq_ids} | parallel -j 1 "samtools view -h \
${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotagged.bam {1} \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp/${barcode_id}_{1}.bam \
&& samtools index ${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp/${barcode_id}_{1}.bam"

# 2. Cluster reads per gene haplotype
cat ${seq_ids} | parallel -j 1 "python $(pwd)/scripts/split_reads.py \
-b ${barcode_id} \
-r {1} \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id} \
${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp/${barcode_id}_{1}.haplotags.tsv \
${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp/${barcode_id}_{1}.bam \
${fastq_file}"

#-------- 5. Assemble consensus from read clusters (spoa) and polish consensus (flye) ----------
## Set flag for Flye's read-based polishing depending on flow cell chemistry
if [[ "${model_config}" == "R9.4.1" ]]; then
	# Use --nano-raw for read-based Flye polishing of R9.4.1 ONT LRAS data
    read_qual=$(echo "--nano-raw")
elif [[ "${model_config}" == "R10.4.1" ]]; then
	# Use --nano-hq for read-based Flye polishing of R10.4.1 ONT LRAS data
    read_qual=$(echo "--nano-hq")
else
    echo "Error: Invalid model flag. Please use -m R9.4.1 for R9.4.1 data or -m R10.4.1 for R10.4.1 LRAS data."
    exit 1
fi

# Check if haplotagged bamfile has read identifiers in the query column before generating draft consensus
parallel '
count=$(samtools view '${hapasm_dir}'/vcf/'${barcode_id}'/'${barcode_id}'_haplotagged.bam {1} | cut -f 1 | wc -l);
if ! [[ $count -ge 2 ]]; then
    region={1}
    echo "No reads for this $region"
else
    for input_fq in $(ls -1 '${hapasm_dir}'/per_gene/{1}/'${barcode_id}'/'${barcode_id}'.{1}.{2}.fastq.gz)
    do
        region={1}
        hap_num={2}
        barcode_id='${barcode_id}'

        base_name=$(basename $input_fq .fastq.gz)
        # hap_num=$(echo "${base_name#$barcode_id.$region.}")
        asm_dir=$(echo "'$hapasm_dir'/per_gene/$region/$barcode_id/assm")
        mkdir -p ${asm_dir}

        spoa -s -l 2 $input_fq > ${asm_dir}/${base_name}.con.fasta

        flye --polish-target ${asm_dir}/${base_name}.con.fasta '${read_qual}' ${input_fq} --iterations 5 --out-dir ${asm_dir}/${hap_num}_flye_polish

        cat ${asm_dir}/${hap_num}_flye_polish/polished_5.fasta | sed -r -e "s@Consensus@${base_name}@g" > ${asm_dir}/${base_name}.fasta

        rm -rf ${asm_dir}/${base_name}.con.fasta*

        rm -rf ${asm_dir}/${hap_num}_flye_polish

        rm -rf ${input_fq}
    done
fi
' :::: $seq_ids ::: h1 h2

#-------------- 6. Directory clean up --------------------

## Remove haplotagged bam of each gene
## rm -rf ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotagged.bam
## rm -rf ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotagged.bam.bai

## Remove TSV file with read ids
rm -rf ${hapasm_dir}/vcf/${barcode_id}/${barcode_id}_haplotags.tsv

## Remove the sample tmp folder
cat $seq_ids | parallel "rm -rf ${hapasm_dir}/per_gene/{1}/${barcode_id}/tmp"

## Remove clair3 log
rm -rf ${hapasm_dir}/vcf/${barcode_id}/log

#------- COMPLETION OF PROTOCOL ------------
echo -e "${barcode_id} haplotypes assembled!!!"
