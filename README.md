# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.

>**R9.4.1 versus R10.4.1 LRAS data:**
>
>By default, this workflow (and accompanying script) uses the `r941_prom_sup_g5014 model` provided with Clair3 for R9.4.1 datasets.
>For users working with R10.4.1 LRAS data, you are encouraged to set the model path for Clair3 in the variant calling step to:
>```bash
>model_path=$(echo "$CONDA_PREFIX/bin/models/r1041_e82_400bps_sup_v500")
>```
>when going through the workflow, or modify **line 72** of the `run_hapasmbl.sh` script to use the `r1041_e82_400bps_sup_v500` model instead.

## Installation of packages using `miniforge`
* Instructions for installing `miniforge` can be found on Jan Kirenz's codelabs (https://kirenz.github.io/codelabs/codelabs/miniforge-setup/#0).
* Packages required to run the HapAsmbl protocol can be installed as shown below.
```bash
# Create environment
mamba create --name hapasm

# Activate environment
mamba activate hapasm

# Install packages
mamba install minimap2 clair3 spoa racon pysam pandas parallel seqkit seqtk bbmap bcftools flye -y

```
## Pipeline: step-by-step
1. Set input variables
```bash
reference="ref.fasta"
fastq_file="F10_VRN2a_FT3.fq.gz"
sample_id="F10"
threads=4
```
2. Read mapping and removal of concatemers - `minimap2`, `samtools` and, `bbmap reformat.sh`
```bash
# Map with Minimap2
minimap2 --MD -a -x map-ont ${reference} ${fastq_file} | samtools sort > ${sample_id}_tmp.bam

# Remove unmapped reads and concatemers
## clipfilter=10 discards reads with more than 10 soft-clipped bases
samtools view -h -F 2308 ${sample_id}_tmp.bam \
| reformat.sh clipfilter=10 in=stdin.bam out=stdout.bam \
| samtools sort > ${sample_id}_clip.bam

# Index bamfile
samtools index ${sample_id}_clip.bam

# Remove tmp bamfile
rm ${sample_id}_tmp.bam
```

3. Variant calling - `clair3`
```bash
# Create working directory for vcf files
mkdir -p vcf

# Set clair3 parameters
platform="ont"
model_path=$(echo "$CONDA_PREFIX/bin/models/r941_prom_sup_g5014")

# run clair3
## options --var_pct_full=1 and --ref_pct_full=1 are recommended for amplicon sequence data
run_clair3.sh \
--bam_fn=${sample_id}_clip.bam \
--ref_fn=${reference} \
--threads=${threads} \
--platform=${platform} \
--model_path=${model_path} \
--output=./vcf/${sample_id} \
--include_all_ctgs \
--sample_name=${sample_id} \
--chunk_size=25000 \
--var_pct_full=1 \
--ref_pct_full=1 \
--no_phasing_for_fa \
--use_whatshap_for_final_output_phasing \
--remove_intermediate_dir

```
4. Read-based phasing of genetic variants into haplotypes - `whatshap`
```bash
# Phase variants
whatshap phase \
-o ./vcf/${sample_id}/${sample_id}_phased.vcf.gz \
--reference ${reference} \
--tag HP \
./vcf/${sample_id}/merge_output.vcf.gz \
${sample_id}_clip.bam \
--indels \
--sample ${sample_id} \
--ignore-read-groups \
--internal-downsampling 23 \
--distrust-genotypes

# Index phased VCF
tabix -f -p vcf ./vcf/${sample_id}/${sample_id}_phased.vcf.gz
```

5. Tag reads from each haplotype in alignment file
* A bamfile (`haplotagged.bam`) will be created in which reads belonging to a haplotype are tagged.
* Reads corresponding to an haplotype will be written to `.tsv` file. 
```bash
# Haplotag reads
whatshap haplotag \
-o ./vcf/${sample_id}/${sample_id}_haplotagged.bam \
--reference ${reference} \
--output-haplotag-list ./vcf/${sample_id}/${sample_id}_haplotags.tsv \
--ignore-read-groups \
--sample ${sample_id} \
--skip-missing-contigs \
./vcf/${sample_id}/${sample_id}_phased.vcf.gz \
${sample_id}_clip.bam

# Index haplotagged bamfile
samtools index ./vcf/${sample_id}/${sample_id}_haplotagged.bam

```

6. For a gene of interest (e.g. _VRN2a_), cluster reads from each haplotype
```bash
target="VRN2a"

# make directory for each gene with sub-directory for each sample
mkdir -p ./${target}/${sample_id}/tmp

# Create an awk expression to extract reads from VRN2a from the .tsv file from previous step
my_awk=$(echo 'BEGIN {OFS = "\t"} /^#/ {print} !/^#/ && $4 ~ var {print}')

cat ./vcf/${sample_id}/${sample_id}_haplotags.tsv \
| awk -v var=VRN2a '$my_awk' > ./${target}/${sample_id}/tmp/${sample_id}_${target}.haplotags.tsv" 

```
7. Read Splitting

   Use TSV file containing reads-haplotype information to split reads according to haplotype

    **Case 1 - Sample is homozygous for reference allele (GT=0/0):**

    * The *${sample_id}_VRN2a.haplotags.tsv* file will be empty, so it is assumed this sample is **HOM_REF**. 

    * In this case, first get alignments of only VRN2a in the haplotagged bam.
    * Then get the read ids from this bam file as shown below:

     ```bash
     # Specify target region
     target="VRN2a"
     
     # Extract reads from only VRN2a in haplotagged bamfile (Very important)
     samtools view -h -o ./${target}/${sample_id}/tmp/${sample_id}_${target}.bam \
     ./vcf/${sample_id}/${sample_id}_haplotagged.bam ${target}

     # Retrieve ids of reads in column 1
     samtools view -h ./${target}/${sample_id}/tmp/${sample_id}_${target}.bam \
     | cut -f 1 > ${sample_id}_${target}_read_ids.txt
     
     ```

   * Afterwards, the IDs are used to exract reads from the sample fastq files and saved into 2 different files `.h1.fastq.gz` and `.h2.fastq.gz` to reflect the original diploid genotype of the sample.

     ```bash
     seqkit grep --pattern-file ${sample_id}_${target}_read_ids.txt ${fastq_file} -o ${sample_id}_${target}.h1.fastq.gz
     seqkit grep --pattern-file ${sample_id}_${target}_read_ids.txt ${fastq_file} -o ${sample_id}_${target}.h2.fastq.gz
     ```

   **Case 2 - Sample is homozygous for alternate allele (GT=1/1):**

   * Although the file isn't empty in this case, the halotype column contains ` none` values, indicating that this sample is **HOM_ALT**.
   * If this is the case, the same steps as in **Case 1** are repeated to retrieve reads from homologous chromosomes.

   **Case 3 - Sample is heterozygous for alternate allele (GT=0/1):**

   * The halotype column in the *${sample_id}_VRN2a.haplotags.tsv* file contains either H1 or H2, which represent the haplotype information of each read. 

   * In this case, the  `whatshap split` program is used to get the read haplotypes from the sample fastq file as shown below.

     ```bash
     # Example syntax: whatshap split --output-h1 h1.fastq.gz --output-h2 h2.fastq.gz reads.fastq.gz haplotypes.txt
     whatshap split \
     --output-h1 ${sample_id}_${target}.h1.fastq.gz \
     --output-h2 ${sample_id}_${target}.h2.fastq.gz \
     ${fastq_file} \
     ./${target}/${sample_id}/tmp/${sample_id}_VRN2a.haplotags.tsv
     ```

An in-house `split_reads.py` python script was used to automate this step. The script was run as shown below.

```bash
mkdir -p cluster_reads

target="VRN2a"

# -r specifies the target region or name of flowering gene
# -b specifies barcode identifier or sample name
# -o specifies output directory to put reads from a haplotype

python ./scripts/split_reads.py \
-b ${sample_id} \
-r ${target} \
-o ./cluster_reads \
./${target}/${sample_id}/tmp/${sample_id}_${target}.haplotags.tsv \
./${target}/${sample_id}/tmp/${sample_id}_${target}.bam \
${fastq_file}

```
8. Generate consensus of reads from each haplotype cluster - `SPOA`
```bash
# Assembly options:
# -s,  --strand-ambiguous: for each sequence pick the strand with the better alignment
# -l, --algorithm:  2 - semi-global

# Example syntax: spoa --strand-ambiguous --algorithm 2 reads.fastq > out.fasta

spoa --strand-ambiguous --algorithm 2 cluster_reads/${sample_id}_${target}.h1.fastq.gz > ${sample_id}_${target}.h1.con.fasta
spoa --strand-ambiguous --algorithm 2 cluster_reads/${sample_id}_${target}.h2.fastq.gz > ${sample_id}_${target}.h2.con.fasta
```

9. Polish consensus with reads - `flye`
```bash
flye --polish-target ${sample_id}_${target}.h1.con.fasta --nano-raw cluster_reads/${sample_id}_${target}.h1.fastq.gz --iterations 5 --out-dir ./
flye --polish-target ${sample_id}_${target}.h2.con.fasta --nano-raw cluster_reads/${sample_id}_${target}.h2.fastq.gz --iterations 5 --out-dir ./
```
10. Head and tail cropping to trim off "foreign" nucleotides (OPTIONAL) - `seqtk`

- For example, to trim off 10 bases from both 5' and 3' ends of contigs,
```bash
seqtk trimfq -b 10 -e 10 ${sample_id}_${target}.h1.fasta > ${sample_id}_${target}.h1.trim.fasta
seqtk trimfq -b 10 -e 10 ${sample_id}_${target}.h2.fasta > ${sample_id}_${target}.h2.trim.fasta
```

#### (OPTIONAL) Using script and sample files provided
To run with one example data (e.g. F10):
```bash
bash ./run_hapasmbl.sh -r ./reference/ref.fasta -f ./fastqs/F10_VRN2a_FT3.fq.gz -o results_test -b F10 -t 8
```
- Final assemblies of VRN2a are located in `results_test/per_gene/VRN2a/F10/assm/F10.VRN2a.h1.fasta` and `results_test/per_gene/VRN2a/F10/assm/F10.VRN2a.h2.fasta`.

To run on multiple samples (e.g. 2 samples):
```bash
parallel -j 1 echo "{} >> sample_ids.txt" ::: F10 F12

cat ./sample_ids.txt | parallel -j 1 "bash ./run_hapasmbl.sh -r ./ref/${reference} -f ./fastqs/{1}_VRN2a_FT3.fq.gz -o results_test -b {1} -t 8"
```


