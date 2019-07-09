#related samples downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/related_samples_vcf/
#unrelated samples downloaded from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

path=/path/to/your/folder


### Family ID, family name, number of individuals (16 for type 1, 2 and 3, and 12 for type 4), sample names

FID=A
family=family$FID
n=`seq 1 16`
samples=`for i in $n; do echo $FID$i; done`
mkdir $path/$family
cd $path/$family

### Extract founders/individuals from a list with vcftools v0.1.14
# Options: --vcf input VCF file, --keep: list of individuals to keep, --out: output file, --recode: generate new file

input=ALL.unrelated_samples.vcf.gz
/soft/bio/vcftools-0.1.14/bin/vcftools --vcfgz $path/$input --keep ../keep_list.txt --out $family.founders --recode


### Remove indels with vcftools v0.1.14
# Options: --vcf input VCF file, --remove-indels: remove indels, --recode: generate new file, --out: output file

/soft/bio/vcftools-0.1.14/bin/vcftools --vcf $family.founders.recode.vcf --remove-indels --recode --out $fmaily.founders.noindel


### Preprare the file for family simulation (remove header, change "|" for "/" and sort by coordinate)

grep ^## $family.founders.noindel.recode.vcf > header1.vcf
grep -v ^# $family.founders.noindel.recode.vcf | sed 's/|/\//g' > $family.founders.txt
sort -k1V -k2n $family.founders.txt > founders_sort.txt 



### Family simulation in R (see file family_simulation.R) run in $path/$family
#########################################################


### Prepare file (set heterozygous as "0/1" and include header with sample ID)

sed 's/1\/0/0\/1/g' family.noheader.vcf > $family.noheader.v2.vcf
echo \#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT $samples | sed 's/ /\t/g' > header2.vcf
cat header2.vcf $family.noheader.v2.vcf > $family.header2.vcf
cat header1.vcf $family.header2.vcf > $family.vcf

rm $family.noheader.v2.vcf
rm header2.vcf
rm $family.header2.vcf

