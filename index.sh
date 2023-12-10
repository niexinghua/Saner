#!/bin/bash
if [ $# -lt 2 ]; then
       echo  "USAGE:sh index.sh ref.fa gff [gtf]"
       exit 0
fi
fa=$1
gff=$2
gtf=$3
echo REF: $fa
echo GFF: $gff
echo GTF: $gtf
scriptdir="$(dirname $0)"

mv $gff $gff.bak
cat $gff.bak|sed  's#gene:##' |sed  's#transcript:##' > $gff



echo -e "\nbuild Index  start:"

echo "RNN CMD:samtools faidx $fa"
samtools faidx $fa

dict=${fa%.*}
echo "RNN CMD: picard CreateSequenceDictionary R=$fa O=$dict.dict"
picard CreateSequenceDictionary R=$fa O=$dict.dict


l=`awk 'BEGIN{a=0}{a=a+$2}END{print a}' $fa.fai`

if [ $l -gt 2000000000 ];then
	echo "RNN CMD: bwa index -a bwtsw $fa"
	bwa index -a bwtsw $fa
else
	echo "RNN CMD: bwa index $fa"
	bwa index $fa

fi

gtfi=""
if [ $gtf ];then
        gtfi=${gtf%.*}
else
	echo -e "\ngtf file not provide, try get gtf from gff:"
        gtfi=${gff%.*}
        echo "RUN CMD: gffread  $gff -T -o "$gtfi.gtf""
        gffread  $gff -T -o "$gtfi.gtf"
		gtf=$gtfi.gtf
fi

echo "build ANNOVAR index"
echo "RUN CMD: gtfToGenePred -genePredExt $gtf unknown_refGene.txt"
gtfToGenePred -genePredExt $gtf unknown_refGene.txt

echo "RUN CMD: retrieve_seq_from_fasta.pl --format refGene --seqfile $fa  unknown_refGene.txt --out unknown_refGeneMrna.fa"
retrieve_seq_from_fasta.pl --format refGene --seqfile $fa --outfile unknown_refGeneMrna.fa unknown_refGene.txt 

