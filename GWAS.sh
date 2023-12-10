#step1 Set Work Path
workdir=/work/my_gwas
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts
tmpdir=$workdir/tmp
export TMPDIR=$tmpdir  
export PATH=$scriptdir:$PATH

GFF=$refdir/chestnut.gff3
REF=$refdir/chestnut.fa
FAI=$refdir/chestnut.fa.fai
GROUP=$datadir/pop_group.txt

#Using emmax for GWAS analysis

plink --vcf $vcf \
     --recode12  --output-missing-genotype 0 \
    --transpose --out snp_var   --set-missing-var-ids @:#  --allow-extra-chr
##Consistency between trait data and SNP genotype samples and convert to Emmax recognition format
perl $scriptdir/sort_trait.pl snp_var.tfam $workdir/05.trait/bio.tsv |awk '{print $1,$2,$NF}' >emmax_bio.tsv
#kinship  ibs-kinship
emmax-kin-intel64 -v  -s -d 10  -o ./pop.kinship.IBS  snp_var
#kinship BN-kinship
emmax-kin-intel64 -v  -d 10  -o ./pop.kinship.BN  snp_var
#pca
awk  '{print $1" "$2" 1 "$3" "$4" "$5" "$6}' ../02.PCA/plink_pca.eigenvec > pca.cov
#Correlation analysis of Emmax model with PCA as covariate
emmax-intel64 -v -d 10 -t snp_var  \
    -p emmax_bio.tsv \
    -k pop.kinship.IBS -c pca.cov  -o emmax.out 
#Organize data and visualize it
awk -F"[:\t]" 'BEGIN{print "CHR\tPOS\tP"}{print $1"\t"$2"\t"$NF}' emmax.out.ps > emmax.pvalue
Rscript $scriptdir/gwas_manhattan_plot.r -i emmax.pvalue -F $FAI -T bio -n bio_emmax_manhattan -c 1.0e-6
Rscript $scriptdir/qq_plot.r  -i emmax.pvalue -n bio_emmax_qq


#Using the LFMM tool  for GWAS analysisin the LEA software package
#Firstly, perform PCA analysis on the bio1-19 data and take the data from PC1 as input data,Then convert the VCF file into an LFMM file
/data/nxh0719/LFMM_SIFT/LFMM_CL_v1.5/bin/LFMM -x /data/nxh0719/LFMM_SIFT/gt.CA4.wild.lfmm -v /data/nxh0719/LFMM_SIFT/pc1.env -K 4 -d 1 -p 20

#Obtain the associated area based on the P-value
ln -s ../07.gwas_tassel/lfmm_pvalue.txt  pvalue.txt
#Filter associated SNPs based on threshold and generate bed files for sorting
awk 'BEGIN{OFS="\t"}NR>1 && $3<1e-7{print $1,$2,$2,$3}' pvalue.txt |sort -k1,1 -k2,2n >pvalue.bed

#Generate a 10k area near the significance P-value
bedtools flank  -i pvalue.bed  -b 10000 -g $FAI  |sort -k1,1 -k2,2n>region.bed
#Merge regions
bedtools merge -i region.bed >region_merged.bed
#Screening genes within associated regions
bedtools  intersect -wo -nonamecheck -F 0.1  -a  region_merged.bed -b gene.gff | awk -F "\t"  '{print $4"\t"$7"\t"$8"\t"$10"\t"$12}' | sort -u  > gene.txt

#Regional LD heat map drawing
LDBlockShow  -InVCF $vcf -OutPut LD -Region Chr1:36000000-37000000 \
    -OutPng -SeleVar 2   
awk '$1=="Chr1" && $2>36000000 && $2<37000000{print}' pvalue.txt >ld.pvalue.txt
LDBlockShow  -Cutline  7 -InVCF $vcf -OutPut LD.p \
    -Region Chr1:36000000-37000000 -OutPng -SeleVar 2 -InGWAS ld.pvalue.txt