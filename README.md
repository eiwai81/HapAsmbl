# HapAsm (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is reference sequence of the amplicons of interest.

## Summary of steps
1. Generate reads using gene-specific primers
2. Pre-assembly processing steps
    - Re-basecalling (dorado)
    - Read splitting at in-read adapters (porechop and/or duplex-tools)
    - Quality filtering (filtlong)
    - Length filtering (filtong)
3. Read mapping/alignment
    - Minimap2
    - LRA
    - NGMLR
    - Vulcan

NOTE: Any of the long-read mappers above can be used. However, minimap2 was used because it was fast and compatible with Clair3 and structural variant callers.

4. Depth normalization (OPTIONAL)
    - samtools
5. SNPs/indel variant calling
    - Longshot - Haplotype-aware variant caller. Calls onls SNPs.
    - Clair3 - Preferred because it also calls indels.
6. Read-based phasing of genetic variants into haplotypes
    - whatshap
7. Cluster unmapped reads into haplotypes
    - whatshap
    - seqkit
8. Generate consensus of reads from each haplotype cluster
    - SPOA
9. Polish consensus with reads
    - Flye
10. Trim head and tail cropping (OPTIONAL)
    - seqtk

## Protocol
### 2. Pre-assembly processing steps
####  2a. Re-basecalling


```
sdfsdfsdfsdf
```
