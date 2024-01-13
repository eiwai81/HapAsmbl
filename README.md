# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.
## Installation of packages with bioconda
All packages (and versions) used in this workflow are explicitly listed in the _hapasmbl_packages.txt_ file. To install them, please use:
```bash
conda create --name hapasmbl --file hapasmbl_packages.txt
conda activate hapasmbl
```
## Summary of workflow
1. Read mapping and removal of concatemers
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
2. Extract reads originating from a flowering gene e.g. CO (_CONSTANS_)
```bash
samtools view -h ${barcode_id}_clip.bam CO -o ${barcode_id}_CO.bam

# Create a sorted and indexed bamfile
samtools sort -o ${barcode_id}_CO.sort.bam ${barcode_id}_CO.bam
samtools index ${barcode_id}_CO.sort.bam

# remove unsorted bamfile
rm ${barcode_id}_CO.bam
```
5. SNPs/indel variant calling
    - `clair3` - Preferred because it also calls indels.
6. Read-based phasing of genetic variants into haplotypes
    - `whatshap`
7. Cluster unmapped reads into haplotypes
    - `whatshap`
    - `seqkit`
8. Generate consensus of reads from each haplotype cluster
    - `SPOA`
9. Polish consensus with reads
    - `flye`
10. Trim head and tail cropping (OPTIONAL)
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


