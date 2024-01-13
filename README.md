# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.
## Installation of packages with bioconda
All packages (and versions) used in this workflow are explicitly listed in the _hapasmbl_packages.txt_ file. To install them, please use:
```bash
conda create --name hapasmbl --file hapasmbl_packages.txt
conda activate hapasmbl
```
## Summary of workflow
1. Set input variables
```bash
reference="ref.fasta"
fastq_file="bc01.fastq"
barcode_id="bc01"
threads=4
```
2. Read mapping and removal of concatemers - `minimap2`, `samtools` and, `bbmap reformat.sh`
```bash
# Map with Minimap2
minimap2 --MD -a -x map-ont ${reference} ${fastq_file} | samtools sort > ${barcode_id}_filt.bam

# Remove unmapped reads and concatemers
## clipfilter=10 discards reads with more than 10 soft-clipped bases
samtools view -h -F 2308 ${barcode_id}_filt.bam \
| reformat.sh clipfilter=10 in=stdin.bam out=stdout.bam \
| samtools sort > ${barcode_id}_clip.bam

# Index bamfile
samtools index ${barcode_id}_clip.bam
```
3. Extract reads originating from a flowering gene e.g. CO (_CONSTANS_) - `samtools`
```bash
# Extract alignments at CO
samtools view -h ${barcode_id}_clip.bam CO -o ${barcode_id}_CO.bam

# Create a sorted and indexed bamfile
samtools sort -o ${barcode_id}_CO.sort.bam ${barcode_id}_CO.bam
samtools index ${barcode_id}_CO.sort.bam

# remove unsorted bamfile
rm ${barcode_id}_CO.bam
```
4. Variant calling - `clair3`
```bash
# Set clair3 parameters
platform="ont"
model_path=$(echo "$CONDA_PREFIX/bin/models/r941_prom_sup_g5014")

# run clair3
run_clair3.sh \
--bam_fn=${barcode_id}_CO.sort.bam \
--ref_fn=${reference} \
--threads=${threads} \
--platform=${platform} \
--model_path=${model_path} \
--output=${barcode_id}_vcf \
--include_all_ctgs \
--sample_name=${barcode_id} \
--bed_fn=CO_region.bed \
--chunk_size=25000 \
--var_pct_full=1 \
--ref_pct_full=1 \
--print_ref_calls \
--snp_min_af=0.01 \
--no_phasing_for_fa \
--use_whatshap_for_final_output_phasing
```
5. Read-based phasing of genetic variants into haplotypes - `whatshap`
```bash
# Phase variants
whatshap phase \
-o ${barcode_id}_CO.phased.vcf.gz \
--reference ${reference} \
--tag HP \
${barcode_id}_vcf/merge_output.vcf.gz \
${barcode_id}_CO.sort.bam \
--indels \
--sample ${barcode_id} \
--ignore-read-groups \
--internal-downsampling 23 \
--distrust-genotypes

# Index phased VCF
tabix -f -p vcf ${barcode_id}_CO.phased.vcf.gz
```
6. Tag reads from each haplotype in alignment file
```bash
# Haplotag reads
whatshap haplotag \
-o ${barcode_id}_CO.haplotagged.bam \
--reference ${reference} \
--output-haplotag-list ${barcode_id}_CO.haplotag_list.tsv.gz \
--ignore-read-groups \
--sample ${barcode_id} \
--skip-missing-contigs \
${barcode_id}_CO.phased.vcf.gz \
${barcode_id}_CO.sort.bam
```
7. Read Splitting - Use TSV file containing reads-haplotype information to split reads according to haplotype

    Case 1 - Sample is homozygous reference allele (GT=0/0):

    * The *.haplotag_list.tsv.gz* file will be empty, so it is assumed this sample is **HOM_REF**. 

    * In this case, the read_IDS are first extracted from the haplotagged bam file produced in the preceding step using the following code:

     ```bash
     # For example, to extract ids from the CO haplotagged bamfile
     samtools view ${barcode_id}_CO.haplotagged.bam | cut -f 1 > ${barcode_id}_CO_read_ids.txt
     
     ```

   * Afterwards, the IDs are used to exract reads from the sample fastq files and saved into 2 different files `.h1.fastq.gz` and `.h2.fastq.gz` to reflect the original diploid genotype of the sample.

     ```bash
     seqkit grep --pattern-file ${barcode_id}_CO_read_ids.txt ${fastq_file} ${barcode_id}_CO.h1.fastq.gz
     seqkit grep --pattern-file ${barcode_id}_CO_read_ids.txt ${fastq_file} ${barcode_id}_CO.h2.fastq.gz
     ```

2. Case 2 - Sample is homozygous for alternate allele (GT=1/1):

   * Although the file isn't empty in this case, the halotype column contains ` none` values which indicates that this sample is **HOM_ALT**. If this is the case, the same steps as in **Case 1** are repeated to retrieve reads from homologous chromosomes.

3. Case 3 - GT=0/1:

   * The halotype column contains either H1 or H2 which represents the haplotype information of each read. 

   * In this case, the  `whatshap split` program to get the read haplotypes as shown below.

     ```bash
     # Example syntax: whatshap split --output-h1 h1.fastq.gz --output-h2 h2.fastq.gz reads.fastq.gz haplotypes.txt
     whatshap split \
     --output-h1 $hapasm_dir/per_gene_bam/VRN1/barcode03/barcode03_VRN1.h1.fastq.gz \
     --output-h2 $hapasm_dir/per_gene_bam/VRN1/barcode03/barcode03_VRN1.h2.fastq.gz \
     $combined_fq/barcode03.trimmed.filt.fastq \
     $hapasm_dir/per_gene_bam/VRN1/barcode03/barcode03_VRN1.haplotag_list.tsv.gz
     
     ```

A crude `in house` python script was used to automate this step. It can be provided on request. The script was run as shown below.

```bash
# For VRN1-VRN3, FT3, FTL9
# split.py was placed in the $working_dir/scripts folder
parallel -j 24 "python ./scripts/split_reads.py \
--region {1} \
--barcode_id {2} \
--haplotagged_bam_dir $hapasm_dir/per_gene_bam \
--fq_dir ./filtered \
$hapasm_dir/per_gene_bam/{1}/{2}/{2}_{1}.haplotag_list.tsv.gz" ::: VRN1 :::: barcode_list.txt

```
10.
11. Generate consensus of reads from each haplotype cluster
    - `SPOA`
12. Polish consensus with reads
    - `flye`
13. Trim head and tail cropping (OPTIONAL)
    - `seqtk`

## Protocol
To run with one example data (e.g. bc01):
```
bash ./run_hapasmbl.sh -r ref/CO.fasta -f ./fastqs/bc01.fastq -o results_test -b bc01 -t 8
```
- Final assemblies are located in `results_test/per_gene/CO/bc01/assm/bc01.CO.h1.fasta` and `results_test/per_gene/CO/bc01/assm/bc01.CO.h2.fasta`.

To run on multiple samples (e.g. 3 samples):
```
parallel -j 1 echo "{} >> sample_ids.txt" ::: bc{01..03}

cat ./sample_ids.txt | parallel -j 1 "bash ./run_hapasmbl.sh -r ref/CO.fasta -f ./fastqs/{1}.fastq -o results_test -b {1} -t 8"
```


