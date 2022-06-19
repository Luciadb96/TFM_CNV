#!/bin/bash

groundtruth_globalset=$1
groundtruth_exome=$2
indet=$3
folder=$4

NCBI_NGS=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/UCSC/NCBI_GET-RM/hom_reg/DeadZone_NGS_exon.bed
NCBI_Sanger=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/UCSC/NCBI_GET-RM/hom_reg/DeadZone_Sanger_exon.bed
GIAB_Wagner=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/UCSC/GIAB/CMRG_GRCh37_exons.bed
blacklistv2=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/UCSC/ENCODE_BLACKLIST/blacklist_v2_hg19.bed
unusual=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/UCSC/UnusualRegionsUCSC/unusual_region_ucsc.bed
pseudogenes=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/PSEUDOGENES/pseudogenes_gencode_v19.bed

repmask=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/REPEATS/repeatMasker.bed
micros=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/REPEATS/microsatellite.bed
segdup=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/REPEATS/genomicSegDups.bed
simprep=/media/lucia/Disco_Duro/ground_truth_NA12878/problematic_regions/REPEATS/simpleRepeats.bed

echo "Ground truth CNVs"
echo "nº of ground truth CNVs (NA12878):  $(wc -l $groundtruth_globalset) "
echo "nº of ground truth CNVs (NA12878) include in an exonic region: $(wc -l $groundtruth_exome)"
echo "nº of ground truth CNVs (NA12878) indeterminate: $(wc -l $indet) "

echo "Problematic Regions"
echo "Nº of highly homologous exonic regions identify in NGS_DeadZone (Mandelker, et al.(2016)): $(wc -l $NCBI_NGS) "
echo "Nº of highly homologous exonic regions identify in Sanger_DeadZone (Mandelker, et al.(2016)): $(wc -l $NCBI_Sanger) "
echo "Nº of exonic regions include in the challenging medically-relevant genes benchmark from Wagner, et al.(2021): $(wc -l $GIAB_Wagner)"
echo "ENCODE blacklist v2: $(wc -l $blacklistv2) "
echo "UCSC Unusual Region: $(wc -l $unusual) "
echo "GENCODE pseudogenes: $(wc -l $pseudogenes) "

echo "Repeat Regions"
echo "Repeating Elements by RepeatMasker: $(wc -l $repmask)"
echo "Microsatellites - Di-nucleotide and Tri-nucleotide Repeats: $(wc -l $micros)"
echo "Segmental duplications: Duplications of >1000 Bases of Non-RepeatMasked Sequence: $(wc -l $segdup)"
echo "Simple Tandem Repeats by TRF: $(wc -l $simprep)" 

mkdir $folder  
#(1) Highly homologous exonic regions identify in NGS_DeadZone
bedtools intersect -wa -a $groundtruth_globalset -b $NCBI_NGS -f 0.3 -c > ./$folder/GTGS_NCBI_NGS.bed
bedtools intersect -wa -a $groundtruth_exome -b $NCBI_NGS -f 0.3 -c > ./$folder/GTE_NCBI_NGS.bed
bedtools intersect -wa -a $indet -b $NCBI_NGS -f 0.3 -c > ./$folder/GTI_NCBI_NGS.bed

#(2) Highly homologous exonic regions identify in Sanger_DeadZone
bedtools intersect -wa -a $groundtruth_globalset -b $NCBI_Sanger -f 0.3 -c > ./$folder/GTGS_NCBI_Sanger.bed
bedtools intersect -wa -a $groundtruth_exome -b $NCBI_Sanger -f 0.3 -c > ./$folder/GTE_NCBI_Sanger.bed
bedtools intersect -wa -a $indet -b $NCBI_Sanger -f 0.3 -c > ./$folder/GTI_NCBI_Sanger.bed

#(3) CMRG benchmark HG002 
bedtools intersect -wa -a $groundtruth_globalset -b $GIAB_Wagner -f 0.3 -c > ./$folder/GTGS_GIAB_Wagner.bed
bedtools intersect -wa -a $groundtruth_exome -b $GIAB_Wagner -f 0.3 -c > ./$folder/GTE_GIAB_Wagner.bed
bedtools intersect -wa -a $indet -b $GIAB_Wagner -f 0.3 -c > ./$folder/GTI_GIAB_Wagner.bed

#(4) ENCODE blacklist v2
bedtools intersect -wa -a $groundtruth_globalset -b $blacklistv2 -f 0.3 -c > ./$folder/GTGS_blacklistv2.bed
bedtools intersect -wa -a $groundtruth_exome -b $blacklistv2 -f 0.3 -c > ./$folder/GTE_blacklistv2.bed
bedtools intersect -wa -a $indet -b $blacklistv2 -f 0.3 -c > ./$folder/GTI_blacklistv2.bed

#(5) UCSC Unusual Region
bedtools intersect -wa -a $groundtruth_globalset -b $unusual -f 0.3 -c > ./$folder/GTGS_unusual.bed
bedtools intersect -wa -a $groundtruth_exome -b $unusual -f 0.3 -c > ./$folder/GTE_unusual.bed
bedtools intersect -wa -a $indet -b $unusual -f 0.3 -c > ./$folder/GTI_unusual.bed

#(6) GENCODE pseudogenes
bedtools intersect -wa -a $groundtruth_globalset -b $pseudogenes -f 0.3 -c > ./$folder/GTGS_pseudogenes.bed
bedtools intersect -wa -a $groundtruth_exome -b $pseudogenes -f 0.3 -c > ./$folder/GTE_pseudogenes.bed
bedtools intersect -wa -a $indet -b $pseudogenes -f 0.3 -c > ./$folder/GTI_pseudogenes.bed

#(7) Repeating Elements by RepeatMasker
bedtools intersect -wa -a $groundtruth_globalset -b $repmask -f 0.3 -c > ./$folder/GTGS_repmask.bed
bedtools intersect -wa -a $groundtruth_exome -b $repmask -f 0.3 -c > ./$folder/GTE_repmask.bed
bedtools intersect -wa -a $indet -b $repmask -f 0.3 -c > ./$folder/GTI_repmask.bed

#(8) Microsatellites
bedtools intersect -wa -a $groundtruth_globalset -b $micros -c > ./$folder/GTGS_micros.bed
bedtools intersect -wa -a $groundtruth_exome -b $micros -c > ./$folder/GTE_micros.bed
bedtools intersect -wa -a $indet -b $micros -c > ./$folder/GTI_micros.bed

#(9) Segmental duplications
bedtools intersect -wa -a $groundtruth_globalset -b $segdup -f 0.3 -c > ./$folder/GTGS_segdup.bed
bedtools intersect -wa -a $groundtruth_exome -b $segdup -f 0.3 -c > ./$folder/GTE_segdup.bed
bedtools intersect -wa -a $indet -b $segdup -f 0.3 -c > ./$folder/GTI_segdup.bed

#(10) Simple Tandem Repeats
bedtools intersect -wa -a $groundtruth_globalset -b $simprep -f 0.3 -c > ./$folder/GTGS_simprep.bed
bedtools intersect -wa -a $groundtruth_exome -b $simprep -f 0.3 -c > ./$folder/GTE_simprep.bed
bedtools intersect -wa -a $indet -b $simprep -f 0.3 -c > ./$folder/GTI_simprep.bed


#bedtools jaccard
mkdir ./$folder/jaccard
#(1) Highly homologous exonic regions identify in NGS_DeadZone
bedtools jaccard -a $groundtruth_globalset -b $NCBI_NGS > ./$folder/jaccard/GTGS_NCBI_NGS.bed
bedtools jaccard -a $groundtruth_exome -b $NCBI_NGS > ./$folder/jaccard/GTE_NCBI_NGS.bed
bedtools jaccard -a $indet -b $NCBI_NGS > ./$folder/jaccard/GTI_NCBI_NGS.bed

#(2) Highly homologous exonic regions identify in Sanger_DeadZone
bedtools jaccard -a $groundtruth_globalset -b $NCBI_Sanger > ./$folder/jaccard/GTGS_NCBI_Sanger.bed
bedtools jaccard -a $groundtruth_exome -b $NCBI_Sanger > ./$folder/jaccard/GTE_NCBI_Sanger.bed
bedtools jaccard -a $indet -b $NCBI_Sanger > ./$folder/jaccard/GTI_NCBI_Sanger.bed

#(3) CMRG benchmark HG002 
bedtools jaccard -a $groundtruth_globalset -b $GIAB_Wagner > ./$folder/jaccard/GTGS_GIAB_Wagner.bed
bedtools jaccard -a $groundtruth_exome -b $GIAB_Wagner > ./$folder/jaccard/GTE_GIAB_Wagner.bed
bedtools jaccard -a $indet -b $GIAB_Wagner > ./$folder/jaccard/GTI_GIAB_Wagner.bed

#(4) ENCODE blacklist v2
bedtools jaccard -a $groundtruth_globalset -b $blacklistv2 > ./$folder/jaccard/GTGS_blacklistv2.bed
bedtools jaccard -a $groundtruth_exome -b $blacklistv2 > ./$folder/jaccard/GTE_blacklistv2.bed
bedtools jaccard -a $indet -b $blacklistv2 > ./$folder/jaccard/GTI_blacklistv2.bed

#(5) UCSC Unusual Region
sort -k1,1 -k2,2n $groundtruth_globalset > global_sort.bed
sort -k1,1 -k2,2n $groundtruth_exome > exome_sort.bed
sort -k1,1 -k2,2n $indet > indet_sort.bed

sort -k1,1 -k2,2n $unusual > unusual_sort.bed
bedtools jaccard -a global_sort.bed -b unusual_sort.bed > ./$folder/jaccard/GTGS_unusual.bed
bedtools jaccard -a exome_sort.bed -b unusual_sort.bed > ./$folder/jaccard/GTE_unusual.bed
bedtools jaccard -a indet_sort.bed -b unusual_sort.bed > ./$folder/jaccard/GTI_unusual.bed

#(6) GENCODE pseudogenes
sort -k1,1 -k2,2n $pseudogenes > pseudogenes_sort.bed
bedtools jaccard -a global_sort.bed -b pseudogenes_sort.bed > ./$folder/jaccard/GTGS_pseudogenes.bed
bedtools jaccard -a exome_sort.bed -b pseudogenes_sort.bed > ./$folder/jaccard/GTE_pseudogenes.bed
bedtools jaccard -a indet_sort.bed -b pseudogenes_sort.bed > ./$folder/jaccard/GTI_pseudogenes.bed

#(7) Repeating Elements by RepeatMasker
sort -k1,1 -k2,2n $repmask > repmask_sort.bed
bedtools jaccard -wa -a global_sort.bed -b repmask_sort.bed > ./$folder/jaccard/GTGS_repmask.bed
bedtools jaccard -wa -a exome_sort.bed -b repmask_sort.bed > ./$folder/jaccard/GTE_repmask.bed
bedtools jaccard -wa -a indet_sort.bed -b repmask_sort.bed > ./$folder/jaccard/GTI_repmask.bed

#(8) Microsatellites
sort -k1,1 -k2,2n $micros > micros_sort.bed
bedtools jaccard -wa -a global_sort.bed -b micros_sort.bed > ./$folder/jaccard/GTGS_micros.bed
bedtools jaccard -wa -a exome_sort.bed -b micros_sort.bed > ./$folder/jaccard/GTE_micros.bed
bedtools jaccard -wa -a indet_sort.bed -b micros_sort.bed > ./$folder/jaccard/GTI_micros.bed

#(9) Segmental duplications
sort -k1,1 -k2,2n $segdup > segdup_sort.bed
bedtools jaccard -wa -a global_sort.bed -b segdup_sort.bed > ./$folder/jaccard/GTGS_segdup.bed
bedtools jaccard -wa -a exome_sort.bed -b segdup_sort.bed > ./$folder/jaccard/GTE_segdup.bed
bedtools jaccard -wa -a indet_sort.bed -b segdup_sort.bed > ./$folder/jaccard/GTI_segdup.bed

#(10) Simple Tandem Repeats
sort -k1,1 -k2,2n $simprep > simprep_sort.bed
bedtools jaccard -wa -a global_sort.bed -b simprep_sort.bed > ./$folder/jaccard/GTGS_simprep.bed
bedtools jaccard -wa -a exome_sort.bed -b simprep_sort.bed > ./$folder/jaccard/GTE_simprep.bed
bedtools jaccard -wa -a indet_sort.bed -b simprep_sort.bed > ./$folder/jaccard/GTI_simprep.bed

rm *_sort.bed

