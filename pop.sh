#phylogenetic tree
#File format conversion
run_pipeline.pl  -Xmx5G -importGuess  $workdir/00.filter/clean.vcf.gz  \
    -ExportPlugin -saveAs supergene.phy -format Phylip_Inter

#Constructing an evolutionary tree using maximum likelihood method
fasttree -nt -gtr  supergene.phy   >  fasttree.nwk

#PCA analysis
#plink
plink --vcf  $workdir/00.filter/clean.vcf.gz --pca 10 --out  plink_pca   \
    --allow-extra-chr --set-missing-var-ids @:#    --vcf-half-call missing
#drawing graph
pca_plink_plot.r -i plink_pca.eigenvec -f $GROUP -g Group --name plink_pca

#STRUCTURE analysis
## filter LD 
plink --vcf  $workdir/00.filter/clean.vcf.gz  --indep-pairwise 50 10 0.2 --out ld   \
    --allow-extra-chr --set-missing-var-ids @:# 
plink --vcf  $workdir/00.filter/clean.vcf.gz  --make-bed --extract ld.prune.in  \
    --out LDfiltered --recode vcf-iid  --keep-allele-order  --allow-extra-chr --set-missing-var-ids @:#  

#Convert to Plink format
vcftools --vcf LDfiltered.vcf --plink \
    --out plink
#Convert to the bed format required by admixture
plink --noweb --file plink  --recode12 --out admixture \
     --allow-extra-chr  --keep-allele-order

#admixture analysis
for k in {2..10};do
    admixture -j2 -C 0.01 --cv admixture.ped $k >admixture.log$k.out
done

#drawing graph 
structure_plot.r  -d ./ -s admixture.nosex 
structure_plot.r  -d ./ -s admixture.nosex -f $GROUP -g Group

#Determine the optimal K
grep "CV error" *out

# LDdecay
for i in cm cs cc ch;do 

    PopLDdecay -InVCF  $workdir/00.filter/clean.vcf.gz \
        -SubPop  ${i}_popid.txt -MaxDist 500 -OutStat ${i}.stat
    Plot_OnePop.pl -inFile ${i}.stat.gz -output ${i}.ld
done
#drawing graph
Plot_MultiPop.pl -inList ld_stat.list -output ld_stat.multi

#Selective sweeping analysis
#Set the input VCF file and calculate the window and step size
gzvcf=$workdir/00.filter/clean.vcf.gz
window=100000
step=10000

#calculate Tajima's D 
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../Cm_popid.txt  --out Cm
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../Cc_popid.txt  --out Cc
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../Cs_popid.txt  --out Cs
vcftools --gzvcf $gzvcf --TajimaD  100000  --keep  ../Ch_popid.txt  --out Ch
#calculate pi
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../NW_popid.txt   --out pi.wild
vcftools  --gzvcf $gzvcf \
    --window-pi $window --window-pi-step  $step  \
    --keep ../SE_popid.txt  --out pi.NW-SE
#calculate fst
vcftools  --gzvcf $gzvcf --fst-window-size $window --fst-window-step $step  \
    --weir-fst-pop  ../NW_popid.txt --weir-fst-pop ../SE_popid.txt --out  Fst.wild.NW-SE
#Fst and Pi joint screening of selected regions
fst_pi_select_sweep.r --fst Fst.wild.cultivated.windowed.weir.fst \
    --pi1 pi.NW.windowed.pi --pi2 pi.SE.windowed.pi --zscore --log2 \
    -A NW -B SE -c 0.05 -n fst-pi.NW-vs-SE  -f pdf
	
#Population dynamics analysis  PSMC
#Chromosome consistency sequences required for generating PSMC
bcftools mpileup -Ou -I -f $REF $workdir/3.map/result/cm-1.sort.dedup.bam| \
bcftools call -c -Ov| vcfutils.pl vcf2fq -d 10 -D 100| gzip - > cm-1.psmc.fq.gz
fq2psmcfa -q20 cm-1.psmc.fq.gz > cm-1.psmc.fa
#Generate chromosome ID list
cat $FAI |grep "Chr"|cut -f1 >chr.list
#Extracting Chromosome Data
 perl /data/nxh0719/scripts_pop1/get_fa_by_id.pl chr.list cm-1.psmc.fa cm-1.Chr.psmc.fa
#psmc analysis
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o cm-1.psmc cm-1.Chr.psmc.fa
##Plot
# -g  The time it takes for this species to reproduce for one generation is - g 30
# -u For the mutation rate of this species
psmc_plot.pl -u 5.2e-8 -g 30  -p cm-1 cm-1.psmc
##bootstrapping
awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length($0)}END{print l}' cm-1.Chr.psmc.fa

splitfa cm-1.Chr.psmc.fa 100000 > cm-1.split.psmcfa
seq 20 | xargs -i echo  psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o round-{}.psmc cm-1.split.psmcfa >boot.sh
parallel -j 10 <boot.sh  # 5个一起并行运行 
cat cm-1.psmc round-*.psmc > cm-1.combined.psmc
