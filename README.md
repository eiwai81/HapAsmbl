# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.
## Installation of packages with bioconda
All packages (and versions) used in this workflow are explicitly listed in the _hapasmbl_packages.txt_ file. To install them or replicate the environment that was used for this project, please use:
```bash
conda create --name hapasmbl --file hapasmbl_packages.txt
conda activate hapasmbl
```
## Pipeline
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

A bed file containing gene region (`CO_region.bed`) flanked by primers was used here.
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
7. Read Splitting

   Use TSV file containing reads-haplotype information to split reads according to haplotype

    **Case 1 - Sample is homozygous reference allele (GT=0/0):**

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

   **Case 2 - Sample is homozygous for alternate allele (GT=1/1):**

   * Although the file isn't empty in this case, the halotype column contains ` none` values which indicates that this sample is **HOM_ALT**. If this is the case, the same steps as in **Case 1** are repeated to retrieve reads from homologous chromosomes.

   **Case 3 - Sample is heterozygous for alternate allele (GT=0/1):**

   * The halotype column in the _.haplotag_list.tsv.gz_ file contains either H1 or H2 which represents the haplotype information of each read. 

   * In this case, the  `whatshap split` program to get the read haplotypes as shown below.

     ```bash
     # Example syntax: whatshap split --output-h1 h1.fastq.gz --output-h2 h2.fastq.gz reads.fastq.gz haplotypes.txt
     whatshap split --output-h1 ${barcode_id}_CO.h1.fastq.gz --output-h2 ${barcode_id}_CO.h2.fastq.gz ${fastq_file} ${barcode_id}_CO.haplotag_list.tsv.gz 
     ```

A crude in-house `split_reads.py` python script was used to automate this step. The script was run as shown below.

```bash
# Specify directory for clustering reads using the option -o

python split_reads.py \
-b ${barcode_id} \
-r CO \
-o path/to/output_dir \
${barcode_id}_CO.haplotag_list.tsv.gz \
${barcode_id}_CO.haplotagged.bam \
${fastq_file}

```
8. Generate consensus of reads from each haplotype cluster - `SPOA`
```bash
# Assembly options:
# -s,  --strand-ambiguous: for each sequence pick the strand with the better alignment
# -l, --algorithm:  2 - semi-global

# Example syntax: spoa --strand-ambiguous --algorithm 2 reads.fastq > out.fasta

spoa --strand-ambiguous --algorithm 2 ${barcode_id}_CO.h1.fastq.gz > ${barcode_id}_CO.h1.con.fasta
spoa --strand-ambiguous --algorithm 2 ${barcode_id}_CO.h2.fastq.gz > ${barcode_id}_CO.h2.con.fasta
```

9. Polish consensus with reads - `flye`
```bash
flye --polish-target ${barcode_id}_CO.h1.con.fasta --nano-raw ${barcode_id}_CO.h1.fastq.gz --iterations 5 --out-dir ./
flye --polish-target ${barcode_id}_CO.h2.con.fasta --nano-raw ${barcode_id}_CO.h2.fastq.gz --iterations 5 --out-dir ./
```
10. Head and tail cropping to trim off "foreign" nucleotides (OPTIONAL) - `seqtk`

#### (OPTIONAL) Using script and sample files provided
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


