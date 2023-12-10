#step1 Set Work Path
workdir=/work/nxh0719
refdir=$workdir/ref
datadir=$workdir/data
scriptdir=$workdir/scripts_req
tmpdir=$workdir/tmp 
REF=$refdir/cm.fasta
BWA_INDEX=$refdir/Cm.fasta
GFF=$refdir/cm.gff3

#Establishing a genomic index
sh $scriptdir/index.sh cm.fasta cm.gff3
#sequence quality control
fastqc $datadir/*.gz  -o $datadir/1.fastqc
#Remove the connector
for line in $(cat chestnut_list.txt)
do
fastp --thread 16 --qualified_quality_phred 10 \
    --unqualified_percent_limit 50 \
    --n_base_limit 10 \
    -i $datadir/${i}_1.fq.gz \
    -I $datadir/${i}_2.fq.gz \
    -o ${i}_1.clean.fq.gz \
    -O ${i}_2.clean.fq.gz \
    --adapter_fasta $workdir/data/illumina_multiplex.fa \
    -h ${i}.html -j ${i}.json
done
#Summary of Quality Control Data Statistics:
python $scriptdir/qc_stat.py -d $workdir/2.data_qc/ -o $workdir/2.data_qc/ -p all_sample_qc


#map bwa 
for line in $(cat chestnut_list.txt)
do
echo  "RUN CMD: bwa mem  $BWA_INDEX  $workdir/work/3.map/bwa/${i}_1.clean.fq.gz \
        $workdir/work/3.map/bwa/${i}_2.clean.fq.gz -t 16 -M \
        -R '@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}' \
        |samtools view -bS -h - > $workdir/3.map/bwa/${i}.bam"

bwa mem  $BWA_INDEX  $workdir/work/3.map/bwa/${i}_1.clean.fq.gz \
        $workdir/work/3.map/bwa/${i}_2.clean.fq.gz -t 64 -M \
        -R "@RG\tID:${i}\tLB:${i}\tPL:ILLUMINA\tSM:${i}" \
        |samtools view -bS -h - > $workdir/work/3.map/bwa/${i}.bam
done

#bam data sorting
for line in $(cat chestnut_list.txt)
do
echo "RUN CMD: picard SortSam -Xmx1g VALIDATION_STRINGENCY=LENIENT I=$workdir/3.map/bwa/${i}.bam \
    O=$workdir/3.map/bwa/${i}.sorted.bam SORT_ORDER=coordinate \
    TMP_DIR=$tmpdir"

picard SortSam -Xmx70g VALIDATION_STRINGENCY=LENIENT I=/work/work/3.map/bwa/${i}.bam \
    O=/work/work/3.map/bwa/${i}.sorted.bam SORT_ORDER=coordinate TMP_DIR=$tmpdir
done

#Remove Duplicates
for line in $(cat chestnut_list.txt)
do
echo "RUN CMD: picard MarkDuplicates -Xmx4g MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 VALIDATION_STRINGENCY=LENIENT  \
    INPUT=$workdir/3.map/bwa/${i}.sorted.bam \
    OUTPUT=$workdir/3.map/result/${i}.sorted.dedup.bam \
    METRICS_FILE=$workdir/3.map/result/${i}.sorted.dedup.metrics"

picard MarkDuplicates -Xmx40g MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512 VALIDATION_STRINGENCY=LENIENT  \
    INPUT=/work/work/3.map/bwa/${i}.sorted.bam \
    OUTPUT=/work/work/3.map/result/${i}.sorted.dedup.bam \
    METRICS_FILE=/work/work/3.map/result/${i}.sorted.dedup.metrics
done

# bam index
for line in $(cat chestnut_list.txt)
 do
    echo "RUN CMD: samtools index ${i}.sorted.dedup.bam"
    samtools index ${i}.sorted.dedup.bam
done


#GATK call SNP

for line in $(cat chestnut_list.txt)
 do
echo "RUN CMD: gatk --java-options '-Xmx40g' HaplotypeCaller -R $REF   \
  I $workdir/3.map/result/${i}.sorted.dedup.bam \
  -O ${i}.g.vcf.gz --max-alternate-alleles 4  --sample-ploidy 2 \
  -ERC GVCF --tmp-dir $tmpdir"

 gatk --java-options "-Xmx4g" HaplotypeCaller -R $REF   \
   -I $workdir/3.map/result/${i}.sorted.dedup.bam \
   -O ${i}.g.vcf.gz --max-alternate-alleles 4  --sample-ploidy 2 \
   -ERC GVCF --tmp-dir $tmpdir 
done

#Merge GVCF
gatk --java-options "-Xmx40g" CombineGVCFs -R $REF \
    --variant p1.g.vcf.gz \
    --variant p2.g.vcf.gz \
    --variant p3.g.vcf.gz \
    --variant p4.g.vcf.gz \
   -O all.g.vcf.gz --tmp-dir $tmpdir

#Convert gvcf to VCF
gatk --java-options "-Xmx4g" GenotypeGVCFs  -R $REF \
  -V all.g.vcf.gz \
  -O all.raw.vcf.gz --tmp-dir $tmpdir

##Variation site quality control, filtering

gatk --java-options "-Xmx4g" VariantFiltration -R $REF -V $workdir/4.snp_indel/GATK/all.raw.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 "  \
    --cluster-window-size 5 --cluster-size 2 --filter-name my_snp_filter -O all.raw.gatk.vcf.gz

zcat all.raw.gatk.vcf.gz |awk '$0~/#/ || ($7 =="PASS"){print $0}' |gzip - > all.raw.gatked.vcf.gz

#vcfutils.pl
# -w INT    SNP within INT bp around a gap to be filtered [3]
# -W INT    window size for filtering adjacent gaps [10]

vcfutils.pl varFilter -w 5 -W 10 "gzip -dc all.raw.gatked.vcf.gz|" |gzip - >all.varFilter.vcf.gz

#vcftools
#--max-missing Exclude sites on the basis of the proportion of missing data 
#(defined to be between 0 and 1, where 0 allows sites that are completely missing 
#and 1 indicates no missing data allowed).

vcftools --gzvcf all.varFilter.vcf.gz --recode --recode-INFO-all --stdout \
    --maf 0.05  --max-missing 0.7  --minDP 4  --maxDP 1000  \
    --minQ 30 --minGQ 50 --min-alleles 2  --max-alleles 2 |gzip - > all.clean.vcf.gz