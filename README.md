# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.
## Installation of packages with bioconda
All packages (and versions) used in this workflow are explicitly listed in the _hapasmbl_packages.txt_ file. To install them, please use:
```
conda create --name hapasm --file hapasmbl_packages.txt
conda activate hapasm
```
## Summary of workflow
1. Generate reads using gene-specific primers
2. Pre-assembly processing steps
    - Re-basecalling
        - `Guppy` (discontinued)
    - Read splitting at in-read adapters (`porechop` and/or `duplex-tools`)
    - Quality filtering (`filtlong`)
    - Length filtering (`filtlong`)
3. Read mapping/alignment
    - `minimap2` (preferred)

**NOTE**: 

Any of the long-read mappers above can be used. However, `minimap2` was used because it was fast and compatible with `clair3` and structural variant callers.

4. Removal of chimeric reads or concatemers
    - `bbmap` (preferred)
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
- All steps of the protocol are in the markdown file `haplotype_assembly_protocol.md`.
- To run with example data:
```
./run_hapasmbl.sh -r ref/CO.fasta -f barcode106_filt.fastq.gz -o results_test -b barcode106 -t 8
```
- To run on multiple samples (e.g. 96 samples):
```
parallel -j 1 echo "{} >> sample_ids.txt" ::: barcode{01..96}

cat ./sample_ids.txt | parallel -j 1 "bash ./run_hapasmbl.sh -r ref/CO.fasta -f {1}_filt.fastq.gz -o results_test -b {1} -t 8"
```

