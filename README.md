# LPSEOsStudy
Codes/Data for LPS/EOs study

StatAnalysisLPS.R was used for analyzing RT-qPCR data from LPS experimental studies

RNA-seq data was analyzed in pipeline as follows: trimgalore (trim and quality) -> hisat (align) -> samtools (sort) -> htseq (counts)
These codes were used in the following order: trimgalore_rnaseq.pbs -> hisat2align.pbs -> bamsorted.pbs -> htseqcount.pbs
Reference genome for alignment: NCBI Refseq (GRCh38).

Count data from htseqcount.pbs was analyzed using DESeq2 (deseq2_code.R) in order to calculate differentially expressed genes 
