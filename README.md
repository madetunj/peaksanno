========================  PEAKS ANNOTATION  ======================

1) USAGE
Executable anywhere as long as the PATHTO is correctly specified.
	```
	PEAKSANNO-main.sh ["MACS PEAK bed"] ["MACS Summit bed"] ["GTF"] ["CHROM SIZES"]
	```

1) INPUT
The package requires the following files
	- MACS peaks bed file
	- MACS peaks summit file
	- GTF having either gene or transcript (prefered) annotation
	- UCSC provided chrom.sizes file

1) OUTPUT

	1) TSS near the center of peaks in __*centerofpeaks_closest.regions.txt*__
	1) peaks overlapping genes regions in __*peaks_within_genebody.regions.txt*__
	1) peaks overlapping promoters in __*peaks_within_promoter.regions.txt*__
	1) peaks overlapping windows in __*peaks_within_window.regions.txt*__
	1) genes identified in previous overlapping regions and comparison of all regions in __*peaks_annotation_comparison.regions.txt*__

1) DEFINITIONS
	- TSS: TSS (transcription start site)
	- promoters: +/- 1kb of TSS
	- window: 10kb upstream --> 3kb downstream of genebody
	- genebody: 1kb upstream of TSS --> TTS/TES (transcription termination/end site)

1) DEPENDENCIES
	* bedtools
	* python3

=================================================================
