 ## A pipeline to assemble haplotypes/alleles of perennial ryegrass flowering genes from multilocus multiplex amplicon sequence data.

1. Software installation
   * The `parallel` software helps to run multiple commands in parallel using shell expansion which allows for parallel processing of several files.

```
conda create -n hapasm python=3.9 -y
conda activate hapasm
conda install minimap2 longshot spoa racon flye samtools bcftools parallel bbmap seqkit seqtk whatshap -y
```

2. Guppy  (version 5.0.16) with GPU was downloaded from Nanopore website.

3. Recall bases from fast5 using commandline Guppy

```bash
# Set absolute path to fast5_pass files generated from GridION
InFolder=$(echo "/Path/to/fast5_pass")

# Set absolute path to write new fastq files to e.g. guppy_out/sup_fast5_pass
out_dir=$(echo "/Path/to/guppy_out/sup_fast5_pass")

# Create output directory
mkdir -p $out_dir

# Path to configuration file to use
supconf=$(echo "/Path/to/Guppy_5.0.16/bin/ont-guppy/data/dna_r9.4.1_450bps_sup.cfg")

parallel -j 4 \
"guppy_basecaller -i ${InFolder}/{1} \
    -s ${out}/{1} \
    -c ${supconf} \
    -x 'cuda:all'" ::: barcode{01..12}
```

4. Combine all re-basecalled reads per barcode into a single file

```bash
# Optional: Create new directory to hold concatenated fastq file per barcode
# set path to project working directory
working_dir=$(echo "/Path/to/working_directory")
combined_fq=$(echo "$working_dir/combined_fastq")

# make directory
mkdir -p $combined_fq

# Concatenate/combine all fastq files in each barcode directory
parallel -j 4 "cat ${out}/{1}/pass/fastq_runid*.fastq ${out}/{1}/pass/fastq_runid*.fastq > $combined_fq/{1}.fastq" ::: barcode{01..12}
```

4. Split adapters within reads and trim terminal adapters with porechop
   * porechop (v 0.2.4)
   * New/updated  porechop repositories: https://github.com/fmfi-compbio/Porechop, https://github.com/dnieuw/Porechop

```bash
# Trim terminal adapters and split adapters in between reads
# Set number of threads to what is available
parallel -j 4 "porechop \
-i $combined_fq/{1}.fastq \
-o $combined_fq/{1}.trimmed.fastq \
--threads 8" ::: barcode{01..12}
```

5. Filter reads based on expected length and average quality score
   * filtlong (v0.2.1)

```bash
# Filtlong parameters
## --min_length: minimum length threshold
## --max_length: maximum length threshold
## --min_mean_q: minimum mean quality threshold

parallel -j 4 "filtlong \
--min_length 1500 \
--max_length 4500 \
--min_mean_q 93 \
$combined_fq/{1}.trimmed.fastq \
> $combined_fq/{1}.trimmed.filt.fastq" ::: barcode{01..12}
```

##### Module 1

Filter out unmapped reads , non-priimary alignments and chimeric/concatemers resulting from direct amplicon-to-amplicon ligation. Mapping requires a reference sequence which can be any other perenial ryegrass reference genomes.

* OPTIONAL (RECOMMENDED):
  * A custom reference sequence containing only genes of interest (*VRN1*, *VRN2a*, *VRN2b*, *FT3*, *FTL9*) was created to quicken total processing time.
  * Reference file was placed in a separate folder in the working directory.

```bash
# Activate conda environment
conda activate hapasm

# Set working directoy
working_dir=$(echo "/Path/to/working_directory")

# Set path to reference
ref=$(echo "$working_dir/ref/flowering_genes.fasta")

# index reference file
samtools faidx $ref

```

###### 1.1 Set other environment variables

```bash
# Set other environment variables
hapasm_dir=$(echo "$working_dir/haplotype_assembly")

# Path to reads filtered by filtlong
input_fq_dir=$(echo "$working_dir/combined_fastq")

```

###### **1.2. Map reads to ref sequence file**

```bash
# Create list of barcode IDs used
## for large number samples, it was best to create a list of barcode ID used
parallel -j 1 echo "{1} >> $working_dir/barcode_list.txt" ::: barcode{01..12}

# make folders for bam, vcfs and fasta file in $hapasm_dir
mkdir -p $hapasm_dir/{bam,per_gene_bam,haplotypes}

# Align the filtered (unmapped) reads to the ref file
## Launch 4 jobs at a time
parallel -j 4 "minimap2 --MD -a -x map-ont \
$ref $input_fq_dir/{1}.trimmed.filt.fastq \
| samtools sort -@ 4 > $hapasm_dir/bam/{}_filt.bam \
&& samtools index $hapasm_dir/bam/{}_filt.bam" :::: barcode_list.txt

# Filter out umapped reads with -F 2308 flag 
# Include -h flag in samtools view to include header information 
## Exclude chimeric reads with samtools with bbmap's reformat.sh
## clipfilter=10 discards reads with more than 10 soft-clipped bases

parallel -j 4 "samtools view -h -F 2308 $hapasm_dir/bam/{}_filt.bam \
| reformat.sh clipfilter=10 in=stdin.bam out=stdout.bam \
| samtools sort > $hapasm_dir/bam/{}_clip.bam \
&& samtools index $hapasm_dir/bam/{}_clip.bam" :::: barcode_list.txt

```

Use the clip.bam files for further analyses.

###### **1.3. Per sample, split reads originating from each flowering gene**

* For each sample, separate reads from *VRN1*, *VRN2a*, *FT3* etc into separate directories.

```bash
# 1. Print genes (or sequence IDs) in the reference and save to a text file
## seqkit was used to retrieve IDs of sequences in reference
## Alternatively, grep -i "^>" $ref > $working_dir/ref/seq_ids.txt can also be used
seqkit seq --name --only-id $ref > $working_dir/ref/seq_ids.txt

# Path to file containing names of sequences in reference file
seq_ids=$(echo "$working_dir/ref/seq_ids.txt")

# make directory for each gene with sub-directory for each sample
parallel "mkdir -p $hapasm_dir/per_gene_bam/{1}/{2}" :::: $seq_ids :::: barcode_list.txt

## Per sample, extract reads that map to each gene region listed in the 
## seq_ids.txt file
# use -h flag in samtools view to include header
parallel -j 4 "samtools view -h $hapasm_dir/bam/{1}_clip.bam {2} -o $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.bam" :::: barcode_list.txt :::: $seq_ids

#----------------------------------------------------------------------------------

# 2. Sort and index the split bams
parallel -j 4 "samtools sort \
-o $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.sort.bam \
$hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.bam \
&& samtools index $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.sort.bam" :::: barcode_list.txt :::: $seq_ids

#----------------------------------------------------------------------------------

# 3. remove unsorted region bamfile
parallel "rm -rf $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.bam" :::: barcode_list.txt :::: $seq_ids

```

##### Module 2

###### 2.1. SNP calling, phasing and haplotyping

* longshot (v0.4.5)
* Clair3 (v 1.0.4)

```bash
# Set path to input and output directories
working_dir=$(echo "/Path/to/working_directory")

# path to current working directory
hapasm_dir=$(echo "$working_dir/haplotype_assembly_050923")

# UPDATE 08/06/2023 - I decided to call variants per region by including
## --use_whatshap_for_final_output_phasing 
## --use_whatshap_for_final_output_haplotagging 
## --bed_fn option
## suggestion by Luo:
## --min_mq=INT   try lowering this to 1 or even 0
## -—no_phasing_for_fa   try enable this
## —var_pct_full=1.0 : try to increase
## —ref_pct_full=1.0 : try to increase
## --min_coverage=6
## --snp_min_af=0.01

#1 VRN1-VRN3, FT3, FTL9
parallel -j 2 "run_clair3.sh \
--bam_fn=${working_dir}/bam/{2}/{2}_clip.bam \
--ref_fn=${ref} \
--threads=2 \
--platform=${platform} \
--model_path=${model_path} \
--output=${working_dir}/variant_analysis/snvs_per_gene/{1}/{2} \
--include_all_ctgs \
--sample_name={2} \
--bed_fn=${bed_path}/{1}.region.bed \
--gvcf \
--chunk_size=25000 \
--var_pct_full=1 \
--ref_pct_full=1 \
--print_ref_calls \
--snp_min_af=0.01 \
--no_phasing_for_fa \
--use_whatshap_for_final_output_phasing" :::: $seq_ids :::: barcode_list.txt

```

##### Module 3

###### 3.1. Phasing VCFs and haplotagging reads 

* whatshap (v1.7)

```bash
#------------- Using clair variant call files ---------------------#
working_dir=$(echo "/Volumes/archive/macknightlab/sequencing_data/nanopore_ampseq_data_all_0623/fq_all_75perc")
ref=$(echo "$working_dir/ref/flowering_genes_kyuss.fasta")

hapasm_dir=$(echo "${working_dir}/haplotype_assembly")
mkdir -p $hapasm_dir

## Alternatively, grep -i "^>" $ref > $working_dir/ref/seq_ids.txt can also be used
#seqkit seq --name --only-id $ref > $working_dir/ref/seq_ids.txt

# Path to file containing names of sequences in reference file
seq_ids=$(echo "$working_dir/ref/seq_ids.txt")

# make directory for each gene with sub-directory for each sample
parallel "mkdir -p $hapasm_dir/per_gene/{1}/{2}" :::: $seq_ids :::: barcode_list.txt

# a. Phasing VCFs and haplotagging reads 
# whatshap to phase SNPs and indels
parallel -j 4 "whatshap phase \
-o $hapasm_dir/per_gene/{2}/{1}/{1}_{2}.phased.vcf.gz \
--reference ${ref} \
--tag HP \
${working_dir}/variant_analysis/snvs_per_gene/{2}/{1}/merge_output.vcf.gz \
${working_dir}/haplotype_assembly/{2}/{1}/{1}_{2}.sort.bam \
--indels \
--sample {1} \
--ignore-read-groups \
--internal-downsampling 23 \
--distrust-genotypes" :::: barcode_list.txt :::: $seq_ids

# b. Index phased compressed vcf.gz file
parallel "tabix -f -p vcf $hapasm_dir/per_gene/{2}/{1}/{1}_{2}.phased.vcf.gz" :::: barcode_list.txt :::: $seq_ids

```

###### 3.2. Haplotag reads in bamfile using phased VCF and generate list of reads belonging to each haplotype group

```bash
# Haplotag reads in bamfile using phased VCF and
# Generate list of reads belonging to a haplotype group 
parallel -j 1 "whatshap haplotag \
-o $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.haplotagged.bam \
--reference ${ref} \
--output-haplotag-list $hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.haplotag_list.tsv.gz \
--ignore-read-groups \
--sample {1} \
--skip-missing-contigs \
$hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.phased.vcf.gz \
$hapasm_dir/per_gene_bam/{2}/{1}/{1}_{2}.sort.bam" :::: barcode_list.txt ::: VRN1

```

###### 3.3. Read Splitting - Use TSV file containing reads-haplotype information to split reads according to haplotype

NOTES:

* Haplotype information of reads in the fastq files will not be produced in the ` .haplotag_list.tsv.gz `  file if:

  * Sample is homozyous for ref allele: Genotype column in VCF is 0/0

  * Sample is homozygous for alt allele: Genotype column in VCF is 1/1
  * No data: Genotype column in VCF is ./.


1. Case 1 - Sample is homozygous reference allele (GT=0/0):
   
   * The *.haplotag.tsv.file* was empty, So it is assumed this sample is **HOM_REF**. 
   
   * In this case, the read_IDS are first extracted from the haplotagged bam file produced in the preceding step using the following code:
   
     ```bash
     # Path to haplotype assembly directory
     hapasm_dir=$(echo "$working_dir/haplotype_assembly")
     
     # For example, to extract ids from the FT3 haplotagged 
     # bam files for barcode01
     samtools view $hapasm_dir/per_gene_bam/VRN1/barcode01/barcode01_VRN1.haplotagged.bam | cut -f 1 > barcode01_VRN1_read_ids.txt
     
     ```
   
   * Afterwards, the IDs are used to exract reads from the sample fastq files and saved into 2 different files `.h1.fastq.gz` and `.h2.fastq.gz` to reflect the original diploid genotype of the sample.
   
     ```bash
     # Path to reads filtered by filtlong
     input_fq_dir=$(echo "$working_dir/combined_fastq")
     
     # Use ids to extract the reads from the filtered sample fastq file
     # For example, to extract FT3 from barcode01 fastq files, see below
     
     ## NOTE: Save reads as h1.fastq.gz to represent FT3 reads from the first homologous chromosome
     seqkit grep --pattern-file barcode01_VRN1_read_ids.txt \
     $combined_fq/barcode01.trimmed.filt.fastq \
     > $hapasm_dir/per_gene_bam/VRN1/barcode01/barcode01.VRN1.h1.fastq.gz
     
     ## NOTE: Save reads as h2.fastq.gz to represent FT3 reads from the second homologous chromosome
     seqkit grep --pattern-file barcode01_VRN1_read_ids.txt \
     $combined_fq/barcode01.trimmed.filt.fastq \
     > $hapasm_dir/per_gene_bam/VRN1/barcode01/barcode01.VRN1.h2.fastq.gz
     
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

##### Module 4

###### 4.1. Consensus generation followed by polishing

1. Use SPOA to build consensus from each read haplotype

   ```bash
   # Create directory in folder containing haplotype fastqs file for assembly and polishing
   parallel -j 4 "mkdir -p $hapasm_dir/per_gene_bam/{1}/{2}/assm" ::: VRN1 :::: barcode_list.txt
   
   # 1. Build consensus from reads
   # 1a. Build consensus of read from haplotype 1
   # -s,  --strand-ambiguous: for each sequence pick the strand with the better alignment
   # -l, --algorithm:  2 -  Use semi-global alignment
   parallel -j 4 "spoa -s -l 2 \
   $hapasm_dir/per_gene_bam/{1}/{2}/{2}_{1}.h1.fastq.gz \
   > $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h1.con.fasta" ::: VRN1 :::: barcode_list.txt
   
   # 1b. Build consensus of read from haplotype 2
   parallel -j 4 "spoa -s -l 2 \
   $hapasm_dir/per_gene_bam/{1}/{2}/{2}_{1}.h2.fastq.gz \
   > $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h2.con.fasta" ::: VRN1 :::: barcode_list.txt
   
   ```

2. Use `flye-polish` to polish consensus

   ```bash
   # 2. Polish consensus assembled with "Flye assembly polisher"
   # flye --polish-target SEQ_TO_POLISH --nano-raw READS --iterations NUM_ITER --out-dir OUTPUTDIR --threads THREADS
   # I used 5 iterations (--iterations 5)
   
   # 2a. polish consensus of first haplotype with reads belonging to this haplotype
   parallel -j 4 "flye \
   --polish-target $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h1.con.fasta \
   --nano-raw $hapasm_dir/per_gene_bam/{1}/{2}/{2}_{1}.h1.fastq.gz \
   --iterations 5 \
   --out-dir $hapasm_dir/per_gene_bam/{1}/{2}/assm/h1_flye_polish" ::: VRN1 :::: barcode_list.txt
   
   # 2b. polish consensus of second haplotype with reads belonging to this haplotype
   parallel -j 4 "flye \
   --polish-target $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h2.con.fasta \
   --nano-raw $hapasm_dir/per_gene_bam/{1}/{2}/{2}_{1}.h2.fastq.gz \
   --iterations 5 \
   --out-dir $hapasm_dir/per_gene_bam/{1}/{2}/assm/h2_flye_polish" ::: VRN1 :::: barcode_list.txt
   
   ```

3. Processing after assembly and polishing

   * Rename sequence header in flye-polished sequences from consensus to barcode.region.hap#.fasta e.g. `barcode01.VRN1.h1.fasta`

   ```bash
   parallel "cat $hapasm_dir/per_gene_bam/{1}/{2}/assm/h1_flye_polish/polished_5.fasta \
   | sed -r -e "s@Consensus@{2}.{1}@g" \
   > $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h1.fasta" ::: VRN1 :::: barcode_list.txt
   
   parallel "cat $hapasm_dir/per_gene_bam/{1}/{2}/assm/h2_flye_polish/polished_5.fasta \
   | sed -r -e "s@Consensus@{2}.{1}@g" \
   > $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h2.fasta" ::: VRN1 :::: barcode_list.txt
   
   ```

4. OPTIONAL - Some house  cleaning

   ```bash
   # Remove assembly directories
   parallel -j 1 "rm -rf $hapasm_dir/per_gene_bam/{1}/{2}/assm/h1_flye_polish" ::: VRN1 :::: barcode_list.txt
   parallel -j 1 "rm -rf $hapasm_dir/per_gene_bam/{1}/{2}/assm/h2_flye_polish" ::: VRN1 :::: barcode_list.txt
   
   # Remove consensus assemblies generated by SPOA
   parallel -j 1 "rm -rf $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h1.con.fasta" ::: VRN1 :::: barcode_list.txt
   parallel -j 1 "rm -rf $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h2.con.fasta" ::: VRN1 :::: barcode_list.txt
   
   ```

###### 4.2. Trim 10 bases from both 5' and 3' ends and place trimmed sequences into the same folder

```bash
# a. Make directory for all VRN1 sequences
parallel -j 1 "mkdir $hapasm_dir/haplotypes/{1}" ::: VRN1

# b. Trim 10 bases at the start and end of each assembly
## Assembly fom haplotype 1
parallel "seqtk trimfq -b 10 -e 10 $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h1.fasta \
> $hapasm_dir/haplotypes/{1}/{2}_{1}.h1.trim.fasta" ::: VRN1 :::: barcode_list.txt

## Assembly fom haplotype 2
parallel "seqtk trimfq -b 10 -e 10 $hapasm_dir/per_gene_bam/{1}/{2}/assm/{2}_{1}.h2.fasta \
> $hapasm_dir/haplotypes/{1}/{2}_{1}.h2.trim.fasta" ::: VRN1 :::: barcode_list.txt

```
