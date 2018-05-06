# GI map tools
This repository contains scripts and notebooks used in "Mapping the Genetic Landscape of Human Cells", Horlbeck et al., in review. At present, these are simply presented as used to process the data in the manuscript, with notebooks serving as examples of their features and use.

## Tripleseq_fastqgz_to_counts.py
Adapted from [https://github.com/mhorlbeck/ScreenProcessing/fastqgz_to_counts.py], this is a command line script that counts reads and barcodes from the sgRNAs and joint barcodes sequenced by the "triple sequencing" approach. Input files are of the format *.fastq.gz and filenames are expected to be sortable such that the read 1, read 2, and read 3 sequence files for each sample are in consecutive order (e.g. sampleA_r1.fastq.gz, sampleA_r2.fastq.gz, sampleA_r3.fastq.gz, sampleB_r1.fastq.gz, sampleB_r2.fastq.gz, sampleB_r3.fastq.gz). The output are three files for each sample containing read counts for every sgRNA pair as determined by counting sgRNAs, barcodes, or reads where the sgRNAs and barcodes were concordent.

## GImap_analysis.py
Collection of functions for computing sgRNA-level and gene-level GI scores and GI correlations, as well as functions for generating key quality control plots.
