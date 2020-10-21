# Xinming Zhuo, xmzhuo@gmail.com
# main script for doing primary and secondary analysis for SNP typing and Haplotyping
# dependency: S4_file.R_script_for_genotype_calling.R, S5_File.R_script_for_haplotype_calling.R in the same folder
# dependency: bwa, samtools, bedtools, bam-readcount, littler , PEAR, bbmap

#!/bin/sh

#feed the designate R script
#access current path or change it to the directory accordingly
R_path=$(pwd)
echo $R_path
#change the file name accordingly
genotype_R="S4_file.R" #R_script_for_genotype_calling
haplotype_R="S5_File.R" #R_script_for_haplotype_calling

#interface for entering inputs
SEPE=`zenity --entry --title="se or pe" --text="1 for se; 2 for pe" --entry-text "2"`
case $? in
	0) 
		echo "\" $SEPE\"  selected.";;
	1) 
		echo "No number entered" ;;
	-1)
		echo "An unexpected error has occured.";;
esac

READLEN=`zenity --entry --title="Read Length of NGS" --text="Input the read length" --entry-text "76"`
case $? in
	0) 
		echo "\" $READLEN\"  selected.";;
	1) 
		echo "No number entered" ;;
	-1)
		echo "An unexpected error has occured.";;
esac

echo $SEPE "*" $READLEN


MINLENGTH=`zenity --entry --title="seed length limit" --text="Enter min length for alignment:" --entry-text "25"`
case $? in
	0) 
		echo "\" $MINLENGTH\"  selected.";;
	1) 
		echo "No number entered" ;;
	-1)
		echo "An unexpected error has occured.";;
esac

echo $MINLENGTH

#!/bin/sh

FILESr=`zenity --file-selection --title="select a reference file" `
case $? in 
	0) 
		echo "\" $FILESr\"  selected.";;
	1)
		echo "No file selected.";;
	-1)
		echo "An unexpected error has occured.";;
esac


for reference in $FILESr
do
   
echo $reference

done

fbname=$(basename "$FILESr" |cut -d. -f1)

echo "ref" $fbname

#!/bin/sh

FILES=`zenity --file-selection --title="select fastq files" --multiple --separator='  '`
case $? in 
	0) 
		echo "\" $FILES\"  selected.";;
	1)
		echo "No file selected.";;
	-1)
		echo "An unexpected error has occured.";;
esac

##extract directory name
DIR=$(dirname "$FILES")
echo "dir" $DIR

echo "$DIR/$fbname.fa"

#copy reference file to target directory
cp $reference "$DIR/$fbname"


# remove R identifier, _R1 , _R2 >_R
filename_f=$(echo $FILES| sed 's/_R./_R/g' )
#remove duplicate names
filename=$(echo $(printf '%s\n' $filename_f|sort -u ))

#for target in $FILES
for target in $filename

do
   
echo $target

#recreate the R1 and R2 name for seeking file
target_R1=$(echo $target| sed 's/_R./_R1_/g' )
target_R2=$(echo $target| sed 's/_R./_R2_/g' )


#Index reference
echo bwa index
bwa index $reference 

#align with mem pe
echo bwa mem
bwa mem -k $MINLENGTH -t 6 $reference $target_R1 $target_R2 > $target.sam #mim length 50

#index reference
echo samtools faidx
samtools faidx $reference

#convert sam to bam
echo samtools view
samtools view -b -S -o $target.bam $target.sam 

#sort bam
echo samtools sort
samtools sort $target.bam $target.sorted 

#report maping percentage
samtools flagstat $target.sorted.bam

#index bam, generate bai
echo samtools index
samtools index $target.sorted.bam

#bcf file
echo samtools mpileup
samtools mpileup -d 1000000 -g -u -D -f $reference $target.sorted.bam > $target.bcf 

#vcf file
echo bcftools view
bcftools view -v -c -g $target.bcf > $target.vcf 

#require preinstall of bam-readcount
#calculate frequency of Nucleotides at each position
echo bam-readcount $reference $target.sorted.bam
~/bam-readcount/bin/bam-readcount -b 15 -d 10000000 -w 10 -l $reference.SNP.bed -f $reference $target.sorted.bam > $target.report.txt

#replace ':' with space
cat $target.report.txt > $target.report.csv
sed -i 's/:/	/g' $target.report.csv

#remapping with pear merge, generated assembled and unassembled fastq. p value 0.001, overlap 50 bp, max aseembled length 250, min assembled length 90
~/pear/bin/pear -y 1G -j 2 -p 0.0001 -v 50 -m 250 -n 90 -f $target_R1 -r $target_R2 -o $target

#reverse compliment the reverse complimented read2
~/bbmap/reformat.sh in=$target.unassembled.reverse.fastq out=$target.unassembled.reverse.rc.fastq rcomp=t ow=t

#bbmap to join non overlapping fastq, with 20N padding
~/bbmap/fuse.sh in1=$target.unassembled.forward.fastq in2=$target.unassembled.reverse.rc.fastq out=$target.nonoverlap.fastq fusepairs=t pad=20 -Xmx200m

#align with mem for overlapped
echo bwa mem
bwa mem -k $MINLENGTH -t 6 $reference $target.assembled.fastq > $target.co.sam 

#align with mem for non overlapped , reduce gap open penalty, reduce penalty for mismatch
bwa mem -k $MINLENGTH -t 6 -A 3 -B 3 -O 3 -E 1 $reference $target.nonoverlap.fastq > $target.cn.sam 

#merge overlap and nonoverlap sam file, use header of co.sam
(grep ^@ $target.sam; for f in "$target.co.sam $target.cn.sam"; do grep -v ^@ $f; done) > $target.c.sam

#require R and litter
#transfer readcount data to R and make file in NextGENe format 
#change 1c to 1d to fix single row issue
Rscript --vanilla $R_path/$genotype_R $target.report.csv $reference.SNP.bed

#read sam, subset sam with SNP bedfile,clustering reads according to similarity
Rscript --vanilla $R_path/$haplotype_R $target.c.sam $reference.SNP.bed $SEPE $READLEN 

rm $target.sam
rm -f $target.c.sam $target.co.sam $target.cn.sam
rm $target.bam
rm $target.assembled.fastq
rm -f $target.unassembled.forward.fastq $target.unassembled.reverse.fastq $target.unassembled.reverse.rc.fastq $target.unassembled.forward.fastq $target.discarded.fastq $target.nonoverlap.fastq

done

echo All done 



