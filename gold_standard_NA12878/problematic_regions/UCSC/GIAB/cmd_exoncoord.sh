#!/bin/bash

# Script to obtain exonic regions include in challenging medically-relevant genes (CMRG) benchmark from Wagner, et al. (2021).

GENEFILE=$1
EXONFILE=$2

#Check if the number of arguments is correct
if [ ! $# -eq 2 ] 
then
	echo -e "The number of arguments is wrong">&2
	exit
fi

OUTPUTFILE=$(basename -s .bed $GENEFILE)_exon.bed

#Print 4th column (gene name list) and search each gene in the exon file.
for gene in $( awk '{ print $4 }' $GENEFILE ); do
	grep -E $gene $EXONFILE >> $OUTPUTFILE ;
done

echo "$OUTPUTFILE file created, with the exonic regions include in CMRG benchmark from Wagner et al (2021)."
