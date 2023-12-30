# HapAsm (Haplotype Assembly)
Protocol for assembling alleles from multilocus long-read amplicon sequence data.

## Summary of steps
1. Generate reads using gene-specific primers
2. Pre-assembly processing steps
    - Re-basecalling
    - Read splitting at in-read adapters
    - Quality filtering
    - Length filtering
3. Read mapping/alignment - Minimap2
    - LRA
    - NGMLR
    - Vulcan*

NOTE: Any of the long-read mappers above can be used. However, minimap2 was used because it was fast and compatible with Clair3 and structural variant callers.

4. Depth normalization (OPTIONAL)
    - samtools
5. SNPs/indel variant calling
    - Longshot - Haplotype-aware variant caller. Calls onls SNPs.
    - Clair3 - Preferred because it also calls indels.
6. Structural variant calling (OPTIONAL)
    - SVIM
    - cuteSV (preferred)
    - sniffles2 
