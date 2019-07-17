# Overview

This repository contains the scripts used for kinship detection using RNA-seq data as described in Blay N, Casas E, Galván-Femenía I, Graffelman J, de Cid R, Vavouri T. Assessment of kinship detection using RNA-seq data. bioRxiv https://doi.org/10.1101/546937.

# Kinship Detection

The kinship_detection_pipeline.sh script does all the analysis from FASTQ files to (i) a .genome file with IBD estimates that can be ploted in R as described in Galván-Femenía et al. (2017), and (ii) a pedigree representation of all the detected familial relationships.

RNA-seq data from all individuals must be stored in FASTQ format in [your_folder]/fastq/.

Additionally, other files are required:
- A BED file with all common SNPs stored at [your_folder]/SNPs/hg19_snp150.bed
- A BED file with RepeatMasker regions stored at [your_folder]/SNPs/hg19_repeatmasker.bed
- A BED file with imprinted genes regions stored at [your_folder]/SNPs/hg19_imprinted.bed
- A HISAT2 index stored at [your_folder]/hisat2_index/genome_tran/
- The reference FASTA file for the human genome stored at [your_folder]/Hsapiens.GRCh37_72.fa
- A .sex file with all samples FID + IID + sex (1=male, 2=female) stored at [your_folder]/

# Family Simulation

The family_simulation.R script simulates the genotypes of the offspring from haplotype data of real unrelated individuals (transformed VCF file). The script prepares the input VCF file with the haplotypes (transformed VCF file) for the family simulation, and after the simulation converts the output to a VCF file with genotype data for all the individuals (real and simulated individuals).

The VCF file with haplotype data of unrelated individuals must be stored at [your_folder]/ALL.unrelated_samples.vcf.gz, and a list of the selected individuals to perform the simulation must be stored at [your_folder]/keep_list.txt.

# Required software

bedtools v2.26.0   
HISAT2 v2.1.0   
samtools v1.7   
vcftools v0.1.14   
plink v1.9   
PRIMUS v1.9.0   
R v3.5   


# Required R packages

data.table   
stringr   
