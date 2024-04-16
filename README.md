# PEAKS ANNOTATION

## DESCRIPTION

Genic annotation of peaks including promoters, gene bodies, gene-centric windows, and proximal genes.

**PeaksAnoo** provides the genic annotation of peaks at promoters, gene bodies, and gene-centric windows.
Annotated regions are collated to give a binary overview of proximal genes, and the peak occupancy 
percentages are graphically presented in a bar plot.
This is a part of the major analysis for ChIPSeq or Cut&Run-Seq using our comprehensive pipeline [SEASEQ](https://github.com/stjude/seaseq/).

## USAGE
``` bash
peaksanno.py -p ["PEAK bed"] [-s ["MACS Summit bed"]] -g ["GTF"] -c ["CHROM SIZES"]
```
## INPUT
The script requires the following files:
1. Peaks bed file
2. GTF/GFF/GFF3 having either gene or transcript (preferred) annotation
3. UCSC provided chromosomal sizes file

An optional summit file from MACS can be included

## OUTPUT
PeaksAnno genes the following files:
1. TSS nearest the center of peaks in **centerofpeaks_closest.regions.txt**.
2. Peaks overlapping genes regions in **peaks_within_genebody.regions.txt**.
3. Peaks overlapping promoters in **peaks_within_promoter.regions.txt**.
4. Peaks overlapping windows in **peaks_within_window.regions.txt**.
5. Peaks identified in previous overlapping regions and comparison of all regions in **peaks_compared_regions.peaks.txt**.
6. Genes identified in previous overlapping regions and comparison of all regions in **peaks_compared_regions.genes.txt**.
7. Distribution graphs in **peaks_compared_regions.distribution.pdf**.

## DEFINITIONS
	- TSS: TSS (transcription start site)
	- promoters: 1kb +/- TSS
	- window: 10kb upstream to 3kb downstream of the gene locus
	- genebody: 1kb upstream of TSS to TES (transcription end site)

## DEPENDENCIES
1. bedtools
2. python3

## CITATION
If you are using ***PEAKSANNO***, please cite its parent paper : 
Adetunji, M.O., Abraham, B.J. SEAseq: a portable and cloud-based chromatin occupancy analysis suite. BMC Bioinformatics 23, 77 (2022). https://doi.org/10.1186/s12859-022-04588-z

## CONTACT

Please use the GitHub issues page for reporting any issues/suggestions (recommended). 

Alternatively, you can e-mail Modupe Adetunji <modupeore.adetunji@stjude.org>

