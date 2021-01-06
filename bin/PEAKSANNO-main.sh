#!/bin/bash
#
# PEAKS ANNOTATION package. Detect peaks based on genome regions
# 01/04/2021

if [ $# -lt 4 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["MACS PEAK bed"] ["MACS Summit bed"] ["GTF"] ["CHROM SIZES"]
  echo ""
  exit 1
fi
#================================================================================
# CONSTANT parameters
BEDFILE=$1
SUMMITFILE=$2
GTF=$3
CHROMSIZES=$4
echo
echo "#############################################"
echo "######       PEAKS ANNOTATION v1       ######"
echo "#############################################"
echo
echo "Required input files provided:"
echo "MACS peaks :  $BEDFILE"
echo "MACS summits :  $SUMMITFILE"
echo "GTF :  $GTF"
echo "UCSC chromsizes :  $CHROMSIZES"
echo
# ================================================================================

tempPrefix=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)"_"`date +%s`

#
# 1 - partition regions from gtf
#
PEAKSANNO-gtf_to_genes.py \
-g $GTF \
-c $CHROMSIZES

#
# 2a - annotate peak regions from bed file for promoters, genebody and window
#
annotationArray=( "promoter" "genebody" "window" )
for args in "${annotationArray[@]}"
do
  # intersect similar peaks with gene regions
  intersectBed -b $BEDFILE -a annotation/$args*-region_genes.gff -wb | \
  cut -f1,4,5,9,10,11,12,13,14,15 | \
  awk -F\\t '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' \
  > $args.$tempPrefix

  # annotate regions found with peak ids
  PEAKSANNO-annotate_regions.py \
  -g $args.$tempPrefix \
  -s $SUMMITFILE \
  -o peaks_within_$args.regions.txt
  export filesArray=$(echo -e "$filesArray\t$args.$tempPrefix")
done

#
# 2b - annotate peak regions from bed file for closest transcripts
#
# extract required columns
cut -f1,4,5,9,10 annotation/TSS-region_genes.gff | \
sort -k1,1 -k2,2n > $tempPrefix-sorted.TSS

# extract peak centers
cat $BEDFILE | awk -F\\t '{print $1 "\t" int(($2+$3)/2) "\t" int(($2+$3)/2) "\t" $4 "\t" $5}' \
| sort -k1,1 -k2,2n \
> $tempPrefix-centerofpeaks.bed

#closest gene features
closestBed -t all -a $tempPrefix-centerofpeaks.bed -b $tempPrefix-sorted.TSS > closest.$tempPrefix

#annotate regions found with peak ids
PEAKSANNO-annotate_regions.py \
-g closest.$tempPrefix \
-s $SUMMITFILE \
-o centerofpeaks_closest.regions.txt
export filesArray=$(echo -e "$filesArray\tclosest.$tempPrefix")

#
# 3 - binary master comparison of all the regions files
#
cut -f 9,10 annotation/original_genes.gff > $tempPrefix.genes.names
PEAKSANNO-include_genes.py \
-g $tempPrefix.genes.names \
-o peaks_annotation_comparison.regions.txt \
-f $filesArray

#
# clean up
#
rm -rf $tempPrefix* $filesArray

echo
echo "COMPLETED"
echo
