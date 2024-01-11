# HapAsmbl (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data. This was originally developed for haplotyping flowering time genes in perennial ryegrass but it can be adapted for any diploid species. It is reference-aided which implies that there is a reference sequence for the amplicons of interest.

## Summary of steps
1. Generate reads using gene-specific primers
2. Pre-assembly processing steps
    - Re-basecalling
        - `Guppy` (discontinued)
        - `dorado`
    - Read splitting at in-read adapters (`porechop` and/or `duplex-tools`)
    - Quality filtering (`filtlong`)
    - Length filtering (`filtlong`)
3. Read mapping/alignment
    - `minimap2` (preferred)
    - `LRA`
    - `NGMLR`
    - `vulcan`

**NOTE**: 

Any of the long-read mappers above can be used. However, `minimap2` was used because it was fast and compatible with `clair3` and structural variant callers.

4. Removal of chimeric reads or concatemers
    - `bbmap` (preferred)
    - `samclip`
5. SNPs/indel variant calling
    - `longshot` - Haplotype-aware variant caller. Calls onls SNPs.
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
All steps of the protocol are in the markdown file `haplotype_assembly_protocol.md`.

