#Identification of harmful mutation sites using SIFT

perl /data/nxh0719/work/SIFT_wild/SIFT4G_Create_Genomic_DB-master/make-SIFT-db-all.pl -config /data/nxh0719/work/SIFT_wild/database/chestnutv1_config.txt

#conda activate sift4g
sift4g -t 96 -d /share/softwares/sift4g/scripts_to_build_SIFT_db/test_files/Gossypium_hirsutum/gene-annotation-src/Gossypium_hirsutum.protein.fasta -q /share/softwares/sift4g/scripts_to_build_SIFT_db/test_files/Gossypium_hirsutum/all_prot.fasta --subst /share/softwares/sift4g/scripts_to_build_SIFT_db/test_files/Gossypium_hirsutum/subst --out /share/softwares/sift4g/scripts_to_build_SIFT_db/test_files/Gossypium_hirsutum/sift4_predictions --sub-results

#SIFT annotation
nohup java -Xmx60G -jar /data/nxh0719/LFMM_SIFT/SIFT/SIFT4G_Annotator.jar  -c -i /data/nxh0719/work/SIFT_wild/database/dbSNP/cm_wild.vcf -d /data/nxh0719/LFMM_SIFT/SIFT/database/chestnut_423_v1 -r /data/nxh0719/work/SIFT_wild/result_cm -t &
#Extracting VCF from different species
nohup bcftools view -S /data/nxh0719/vcf/Cm_popid.txt /data/nxh0719/vcf/426wild/426wild.china0.2cq50.snp.vcf.gz -Ov  |gzip - > /data/nxh0719/work/SIFT_wild/database/dbSNP/cc_wild.vcf.gz &

#Statistically identifying the types and number of sites for harmful mutations
perl /data/nxh0719/work/SIFT_wild/vcf_trans2numtype.pl selected_snps.vcf