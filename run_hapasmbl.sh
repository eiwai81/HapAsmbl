#!/usr/bin/env bash

# dead programs tell no lies
set -ueo pipefail

# Get options
# Usage information
# check whether user had supplied -h or --help . If yes display usage
if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
	echo "Usage: $0 [-r reference_file] [-f fastq_file] [-o output_dir] [-b barcode_id] [-t threads]"
	exit 0
fi

# Get commandline arguments
while getopts r:f:o:b:t: flag
do
    case "${flag}" in
        r) reference=${OPTARG};;
        f) fastq_file=${OPTARG};;
        o) out_dir=${OPTARG};;
        b) barcode_id=${OPTARG};;
        t) threads=${OPTARG};;
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

# Create working directory
mkdir -p ${hapasm_dir}/bam/${barcode_id}

# Create sub-directory for each amplicon 
parallel "mkdir -p ${hapasm_dir}/per_gene/{}" :::: ${seq_ids}

######## MODULE 1 ##########################

# ---------- 1.2 Map reads to reference ---------------
echo "Mapping filtered (unmapped) reads to the reference file"
minimap2 --MD -a -x map-ont ${reference} ${fastq_file} | samtools sort > ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_filt.bam

# Sort indexed bam file
samtools index ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_filt.bam

echo "Removing unmapped and chimeric reads"
# Filter out umapped reads with -F 2308 flag 
# Include -h flag in samtools view to include header information 
## Exclude chimeric reads with samtools with bbmap's reformat.sh
## clipfilter=10 discards reads with more than 10 soft-clipped bases

samtools view -h -F 2308 ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_filt.bam \
| reformat.sh clipfilter=10 in=stdin.bam out=stdout.bam \
| samtools sort > ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam

samtools index ${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam

# ----------- 1.3 Per sample, split reads originating from each flowering gene ------
# make directory for each gene with sub-directory for each sample
cat ${seq_ids} | parallel -j 1 "mkdir -p ${hapasm_dir}/per_gene/{1}/${barcode_id}"

# 1. per sample, extract reads that map to each gene region listed in the seq_ids.txt file
# use -h flag in samtools view to include header
cat ${seq_ids} | parallel -j 1 "samtools view -h \
${hapasm_dir}/bam/${barcode_id}/${barcode_id}_clip.bam {1} \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.bam"

#----------------------
# 2. sort and index the split bams
cat ${seq_ids} | parallel -j 1 "samtools sort \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.sort.bam \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.bam \
&& samtools index ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.sort.bam"

#------------------------

# 3. remove unsorted region bamfile
cat ${seq_ids} | parallel "rm -rf ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.bam"

######## MODULE 2 ##########################

# ---------- 2.1 Call variants - clair3  ---------------
# Call variants per region by including with clair3
platform="ont"
model_path=$(echo "$CONDA_PREFIX/bin/models/r941_prom_sup_g5014")

# path to bed
bed_path=$(echo "$(pwd)/ref")

cat ${seq_ids} | parallel -j 1 "run_clair3.sh \
--bam_fn=${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.sort.bam \
--ref_fn=${reference} \
--threads=${threads} \
--platform=${platform} \
--model_path=${model_path} \
--output=${hapasm_dir}/per_gene/{1}/${barcode_id} \
--include_all_ctgs \
--sample_name=${barcode_id} \
--bed_fn=${bed_path}/{1}.region.bed \
--gvcf \
--chunk_size=25000 \
--var_pct_full=1 \
--ref_pct_full=1 \
--print_ref_calls \
--snp_min_af=0.01 \
--no_phasing_for_fa \
--use_whatshap_for_final_output_phasing"

# ---------- 3.1 Phase variants - whatshap ---------------
cat ${seq_ids} | parallel -j 1 "whatshap phase \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.phased.vcf.gz \
--reference ${reference} \
--tag HP \
${hapasm_dir}/per_gene/{1}/${barcode_id}/merge_output.vcf.gz \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.sort.bam \
--indels \
--sample ${barcode_id} \
--ignore-read-groups \
--internal-downsampling 23 \
--distrust-genotypes"

# Index phased compressed vcf.gz file
cat ${seq_ids} | parallel "tabix -f -p vcf ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.phased.vcf.gz"

#------- 3.2 Haplotag reads in bamfile using phased VCF and ----------------------------
#-------     generate list of reads belonging to each haplotype group

cat ${seq_ids} | parallel -j 1 "whatshap haplotag \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.haplotagged.bam \
--reference ${reference} \
--output-haplotag-list ${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.haplotag_list.tsv.gz \
--ignore-read-groups \
--sample ${barcode_id} \
--skip-missing-contigs \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.phased.vcf.gz \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.sort.bam"


#-------- 3.3 Cluster reads from each haplotype ---------

cat ${seq_ids} | parallel -j 1 "python split_reads.py \
-b ${barcode_id} \
-r {1} \
-o ${hapasm_dir}/per_gene/{1}/${barcode_id} \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.haplotag_list.tsv.gz \
${hapasm_dir}/per_gene/{1}/${barcode_id}/${barcode_id}_{1}.haplotagged.bam \
${fastq_file}"

#-------- 4.1 Assemble and polish consensus with spoa and flye ----------

parallel -j 24 "bash hapasm_module3.sh {1} {2} {3} {4} \
${hapasm_dir}/per_gene/{1}/{2}/assm \
${hapasm_dir}/per_gene/{1}/{2}/{2}.{1}.{3}.fastq.gz" :::: $seq_ids ::: ${barcode_id} ::: h1 h2 ::: 2




















