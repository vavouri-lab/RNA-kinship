
path=/path/to/your/folder
#provide a name for the family or group of individuals
family=CU1463
#individuals to be analyzed
samples="SRR1258217
SRR1258218
SRR1258219
SRR1258220
SRR1258221
SRR1258222
SRR1258223
SRR1258224
SRR1258225
SRR1258226
SRR1258227
SRR1258228
SRR1258229
SRR1258230
SRR1258231
SRR1258232
SRR1258233"

####### Selection of SNPs from common hg19 snp150 (UCSC dbSNP) with bedtools intersect v2.26.0
# Options: -v: report those entries in -a that have no overlaps with -b

# not in repeats (from RepeatMasker)

intersectBed -v -a $path/SNPs/hg19_snp150.bed -b $path/SNPs/hg19_repeatmasker.bed > $path/SNPs/hg19_snp150_repeatmasker.bed

# not in imprinted genes (from geneimprint)

intersectBed -v -a $path/SNPs/hg19_snp150_repeatmasker.bed -b $path/SNPs/imprinted_genes.bed > $path/SNPs/hg19_snp150_rmi.bed


####### Kinship detection pipeline (from .fastq)

mkdir $path/$family
### Map to the reference genome with HISAT2
# Options: -p 7: use 7 threads, -q: input files are FASTQ, -x: index(ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch37_tran.tar.gz), -1: first pair, -2: second pair

### Convert from SAM to BAM with samtools view
# Options: -b: output BAM, -S: input format

for ID in $samples
do hisat2 -p 7 -q -x $path/hisat2_index/genome_tran -1 $path/fastq/$ID"_1.fastq" -2 $path/fastq/$ID"_2.fastq" | samtools view -bS - > $path/$family/$ID.bam
done


### Remove duplicates with samtools markdup v1.7
# Options: -@6: use 6 additional threads, sort -n: sort by name, -o: output file, fixmate -m: add MC and mate score tag, sort: sort by coordinates, markdup -r: remove duplicates

for ID in $samples
do
samtools sort -@ 6 -n -o $path/$family/namesort_$ID.bam $path/$family/$ID.bam
samtools fixmate -@ 6 -m $path/$family/namesort_$ID.bam $path/$family/fixmate_$ID.bam
samtools sort -@ 6 -o $path/$family/positionsort_$ID.bam $path/$family/fixmate_$ID.bam
samtools markdup -@ 6 -r $path/$family/positionsort_$ID.bam $path/$family/markdup_$ID.bam
rm $path/$family/namesort_$ID.bam
rm $path/$family/fixmate_$ID.bam
rm $path/$family/positionsort_$ID.bam
done


### Call genotypes of the selection of SNPs with samtools mpileup v1.7 + bcftools call v1.4
# Options: -A: do not discard anomalous read pairs, -q 4: uniquely mapping reads, -t AD,DP: add AD and DP tags, -l: positions file, -f: fasta reference file, -g: BCF format, call -m: multiallelic caller, -O b: output BCF -f GQ: add GQ field

samtools mpileup -A -q 4 -t AD,DP -l $path/SNPs/hg19_snp150_rmi.bed -f /biodata/indices/species/ENSEMBL_nr/Hsapiens/Hsapiens.GRCh37_72.fa $path/$family/markdup_*.bam -g | bcftools call -m - -O b -f GQ > $path/$family/hg19_snp150_markdup_rmi_hisat_$family.q4.bcf


### Filter SNPs (GQ >= 20, DP >= 10) with vcftools v0.1.14
# Options: --bcf: input BCF file, --out: output file, --minGQ 20: exclude genotypes with a Genotype Quality lower than 20, --minDP 10: exclude genotypes with a Depth lower than 10, --recode: generate a new file, --recode-INFO-all: keep all INFO fields

vcftools --bcf $path/$family/hg19_snp150_markdup_rmi_hisat_$family.q4.bcf --out $path/$family/hg19_snp150_markdup_rmi_hisat_$family"_"q4_GQ20_DP10 --minGQ 20 --minDP 10 --recode --recode-INFO-all


### Pairwise comparisons with PLINK2 v1.9 to obtain IBD estimates
# Options: --vcf: input VCF, --maf 0.3: exclude genotypes with a MAF lower than 0.3, --make-bed: create .bed+.bim+.fam files, --out: output file, --bfile: input .bed+.bim+.fam files, --genome full: create a .genome IBD estimates file with all fields

cd $path/$family/
input=hg19_snp150_markdup_rmi_hisat_$family"_"q4_GQ20_DP10.recode.vcf
maf=3
plink --vcf $input --maf 0.$maf --make-bed --out $input.maf$maf
plink --bfile $input.maf$maf --genome full --out $input.maf$maf


### Pedigree representation with PRIMUS v1.9.0
# Options: -p: input plink .genome file, -t 0.2 (or -t 0.375): minimum level of relatedness for related individuals (otherwise they will be considered unrelated), --sex_file: a file with all samples FID + IID + sex (1=male, 2=female), -o: output directory (will be created)

run_PRIMUS.pl -p $path/$family/$input.maf$maf.genome -t 0.2 --sex_file $path/$family/$family.sex -o $path/$family/primus.$family





