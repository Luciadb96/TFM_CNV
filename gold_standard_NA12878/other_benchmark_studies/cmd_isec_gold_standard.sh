#!/bin/bash

GT=$1
exome=$2
indet=$3

paper1=gold_standard_NA12878_Sun.bed
paper2=gold_standard_NA12878_Haraksingh.bed
paper3=gold_standard_NA12878_Zhang.bed
paper4=gold_standard_NA12878_Parikah.bed

folder=$4

echo " nº of CNV intervals in the ground truth:  $(wc -l $GT) "
echo " nº of exon-CNV intervals in the ground truth:  $(wc -l $exome) "
echo " nº of indeterminate CNV intervals:  $(wc -l $indet) "
echo " nº of CNV calls in Sun (NA12878):  $(wc -l $paper1) "
echo " nº of CNV calls in Haraksingh (NA12878):  $(wc -l $paper2) "
echo " nº of CNV calls in Zhang (NA12878):  $(wc -l $paper3)"
echo " nº of CNV calls in Parikah (NA12878):  $(wc -l $paper4)"

mkdir ./$folder
mkdir ./$folder/result/
mkdir ./$folder/jaccard/

#Bedtools intersect
bedtools intersect -wa -a $GT -b $paper1 -f 0.3 -r -c > ./$folder/result/CGT_paper1_isec.bed
bedtools intersect -wa -a $exome -b $paper1 -f 0.3 -r -c > ./$folder/result/Cexome_paper1_isec.bed
bedtools intersect -wa -a $indet -b $paper1 -f 0.3 -r -c > ./$folder/result/Cindet_paper1_isec.bed

bedtools intersect -wa -a $GT -b $paper2 -f 0.3 -r -c > ./$folder/result/CGT_paper2_isec.bed
bedtools intersect -wa -a $exome -b $paper2 -f 0.3 -r -c > ./$folder/result/Cexome_paper2_isec.bed
bedtools intersect -wa -a $indet -b $paper2 -f 0.3 -r -c > ./$folder/result/Cindet_paper2_isec.bed

bedtools intersect -wa -a $GT -b $paper3 -f 0.3 -r -c > ./$folder/result/CGT_paper3_isec.bed
bedtools intersect -wa -a $exome -b $paper3 -f 0.3 -r -c > ./$folder/result/Cexome_paper3_isec.bed
bedtools intersect -wa -a $indet -b $paper3 -f 0.3 -r -c > ./$folder/result/Cindet_paper3_isec.bed

bedtools intersect -wa -a $GT -b $paper4 -f 0.3 -r -c > ./$folder/result/CGT_paper4_isec.bed
bedtools intersect -wa -a $exome -b $paper4 -f 0.3 -r -c > ./$folder/result/Cexome_paper4_isec.bed
bedtools intersect -wa -a $indet -b $paper4 -f 0.3 -r -c > ./$folder/result/Cindet_paper4_isec.bed

#bedtools jaccard
bedtools jaccard -a $GT -b $paper1 > ./$folder/jaccard/GT_Sun_jac.bed
bedtools jaccard -a $exome -b $paper1 > ./$folder/jaccard/exome_Sun_jac.bed
bedtools jaccard -a $indet -b $paper1 > ./$folder/jaccard/indet_Sun_jac.bed

bedtools jaccard -a $GT -b $paper2 > ./$folder/jaccard/GT_Harak_jac.bed
bedtools jaccard -a $exome -b $paper2 > ./$folder/jaccard/exome_Harak_jac.bed
bedtools jaccard -a $indet -b $paper2 > ./$folder/jaccard/indet_Harak_jac.bed

bedtools jaccard -a $GT -b $paper3 > ./$folder/jaccard/GT_Zhang_jac.bed
bedtools jaccard -a $exome -b $paper3 > ./$folder/jaccard/exome_Zhang_jac.bed
bedtools jaccard -a $indet -b $paper3 > ./$folder/jaccard/indet_Zhang_jac.bed

bedtools jaccard -a $GT -b $paper4 > ./$folder/jaccard/GT_Parikah_jac.bed
bedtools jaccard -a $exome -b $paper4 > ./$folder/jaccard/exome_Parikah_jac.bed
bedtools jaccard -a $indet -b $paper4 > ./$folder/jaccard/indet_Parikah_jac.bed
