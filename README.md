#PEAKS ANNOTATION

## DESCRIPTION

Genic annotation of peaks including promoters, gene bodies, gene-centric windows, and proximal genes.

PEAKSANNO generates the genic annotation of peaks at promoters, gene bodies and gene-centric windows.
Annotated regions are collated to provide a binary overview of proximal genes, and the peaks occupancy 
percentages are graphically presented in a bar plot.
This is a part of the major analysis for ChIPSeq or Cut&Run-Seq using our comprehensive pipeline [SEASEQ](https://github.com/stjude/seaseq/).

## USAGE.
	```
	peaksanno.py -p ["MACS PEAK bed"] -s ["MACS Summit bed"] -g ["GTF"] -c ["CHROM SIZES"]
	```
## INPUT
The package requires the following files
	- Peaks bed file
	- GTF/GFF/GFF3 having either gene or transcript (prefered) annotation
	- UCSC provided chrom.sizes file
Optional summit file from MACS can be included

## OUTPUT

	1) TSS nearest the center of peaks in __*centerofpeaks_closest.regions.txt*__
	1) Peaks overlapping genes regions in __*peaks_within_genebody.regions.txt*__
	1) Peaks overlapping promoters in __*peaks_within_promoter.regions.txt*__
	1) Peaks overlapping windows in __*peaks_within_window.regions.txt*__
	1) Peaks identified in previous overlapping regions and comparison of all regions in __*peaks_compared_regions.peaks.txt*__
	1) Genes identified in previous overlapping regions and comparison of all regions in __*peaks_compared_regions.genes.txt*__
	1) Distribution graphs in __*peaks_compared_regions.distribution.pdf*__

## DEFINITIONS
	- TSS: TSS (transcription start site)
	- promoters: 1kb +/- TSS
	- window: 10kb upstream to 3kb downstream of genebody
	- genebody: 1kb upstream of TSS to TES (transcription end site)

## DEPENDENCIES
	* bedtools
	* python3


##Contact

Please use the GitHub issues page for reporting any issues / suggestions (recommended). 

Alternatively, you can e-mail Modupe Adetunji <modupeore.adetunji@stjude.org>

