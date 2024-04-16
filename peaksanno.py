#!/usr/bin/env python3
'''
========================  PEAKS ANNOTATION  ======================

USAGE
All inclusive script for Annotation of ChIP-Seq peaks.
This is included in the SEASEQ pipeline.
    ```
    python PEAKSANNO.py -g <GTF/GFF/GFF3> -c <CHROM SIZES> -p <MACS PEAKS bed> -s <MACS SUMMITS bed>
	```

The script input requirements are :
- GTF : GTF/GFF/GFF3 having either gene or transcript (prefered) annotation
- CHROM SIZES : chrom.sizes file
- MACS PEAKS bed : MACS peaks bed file
- MACS SUMMITS bed : MACS peaks summit file

Output files provided are :
	1) TSS nearest the center of peaks in __*centerofpeaks_closest.regions.txt*__
	2) peaks overlapping genes regions in __*peaks_within_genebody.regions.txt*__
	3) peaks overlapping promoters in __*peaks_within_promoter.regions.txt*__
	4) peaks overlapping windows in __*peaks_within_window.regions.txt*__
	5) peaks identified in previous overlapping regions and comparison of
           all regions in __*peaks_compared_regions.peaks.txt*__
	6) genes identified in previous overlapping regions and comparison of
           all regions in __*peaks_compared_regions.genes.txt*__
	7) distribution graphs in __*peaks_compared_regions.distribution.pdf*__

Definitions of annotation regions:
	- TSS: TSS (transcription start site)
	- promoters: 1kb +/- TSS
	- window: 10kb upstream --> 3kb downstream of genebody
	- genebody: 1kb upstream of TSS --> TES (transcription end site)

Dependencies of script:
	* bedtools
	* python3

==================================================================
'''

import os
import re
import argparse
import subprocess
import matplotlib.pyplot as plt
from pathlib import Path
import collections
from collections import defaultdict

if not os.path.exists('annotation'):
    os.makedirs('annotation')

def parse_genelocations(chromz, results, up_dist, down_dist, tss=True):
    """ Parse genomic regions
    Args:
        chromz (dict) : Chromosomal sizes
        results (string) : Individual gene coordinates
        up_dist (int) : genomic distance upstream / +
        down_dist (int) : genomic distance downstream / -
        tss (boolean) : TSS or entire gene (true if tss)
    """
    lines = results.split("\t")
    lines[3], lines[4] = int(lines[3]), int(lines[4])

    if lines[6] == "+":
        end_region = 4 if tss else 3
        end_position = lines[end_region] + down_dist
        start_position = lines[3] - up_dist
    elif lines[6] == "-":
        end_region = 3 if tss else 4
        start_position = lines[end_region] - down_dist
        end_position = lines[4] + up_dist

    start_position = max(start_position, 1)
    end_position = min(end_position, int(chromz[lines[0]]))

    return "{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]), start_position, end_position, "\t".join(lines[5:]))


def gtf_to_genes(chrom_sizes_file, gtf_file):
    """ GTF to GENES
    Args:
        chrom_sizes_file (file) : Chromosomal sizes file
        gtf_file (file) : GTF file
    """
    #read ucsc chrom sizes, gtf
    chrom_sizes_dict = {line.split('\t')[0]: line.split('\t')[1].rstrip("\n") for line in open(chrom_sizes_file)}
    haschr = any(chrom.startswith('chr') for chrom in chrom_sizes_dict)  #check 'chr' prefix in gtf file

    feature_dict = {line.split("\t")[2]: line.split("\t")[2] for line in open(gtf_file) if not line.startswith('#')}

    feature = "transcript" if 'transcript' in feature_dict else "gene"
    print("NOTE :\tFeature type used is '%s'" % feature)

    #parsing gene annotation file, GFF, GFF3, GTF
    with open(gtf_file, 'r') as gtf_input:
        for gtf_col in gtf_input:
            if not gtf_col.startswith('#'):
                gtf_col.split("\t")[0] = "chr" + gtf_col.split("\t")[0] if haschr and not gtf_col.split("\t")[0].startswith('chr') else gtf_col.split("\t")[0][3:] if not haschr and gtf_col.split("\t")[0].startswith('chr') else gtf_col.split("\t")[0]
                if gtf_col.split("\t")[2] == feature:
                    record = gtf_col.split("\t")[8].split(';')
                    if gtf_file.split('.')[-1] in ['gff', 'gff3']:
                        transcript_id = [s for s in record if "transcript_id" in s][0] if re.search(";transcript_id", gtf_col.split("\t")[8]) else record[0]
                        gene_name = [s for s in record if "gene_name" in s][0] if re.search(";gene_name", gtf_col.split("\t")[8]) else record[0]
                        transcript_id = re.sub('[\"\;]', '', transcript_id.split('=')[-1]) #clean transcript_id
                        gene_name = re.sub('[\"\;]', '', gene_name.split('=')[-1]) #clean gene_name
                    elif gtf_file.split('.')[-1] == 'gtf':
                        transcript_id = [s for s in record if "transcript_id " in s][0] if re.search("; transcript_id", gtf_col.split("\t")[8]) else record[0]
                        gene_name = [s for s in record if "gene_name " in s][0] if re.search("; gene_name", gtf_col.split("\t")[8]) else record[0]
                        transcript_id = re.sub('[\"\;]', '', transcript_id.split(' ')[-1]) #clean transcript_id
                        gene_name = re.sub('[\"\;]', '', gene_name.split(' ')[-1]) #clean gene_name

                    results = "{0}\t{1}\t{2}\t{3}".format(gtf_col.split("\t")[0], "\t".join(gtf_col.split("\t")[1:8]), gene_name, transcript_id)

                    try:
                        if chrom_sizes_dict[gtf_col.split("\t")[0]]:
                            with open('annotation/original_genes.gff', 'w') as original_output:
                                original_output.write(results + "\n")
                                final_output = parse_genelocations(chrom_sizes_dict, results, 1000, 1000, False)
                                for file_name, up, down, tss in [('promoter-region_genes.gff', 10000, 3000, True), ('window-region_genes.gff', 1000, 1000, False), ('TSS-region_genes.gff', 0, 0, False), ('genebody-region_genes.gff', 1000, 0, True)]:
                                    with open(f'annotation/{file_name}', 'w') as file:
                                        file.write(final_output)
                    except KeyError:
                        continue


def annotate_regions(region_file, output_file, macs_summits_file=None):
    """ Annotate Peaks based on genomic features
    Args:
        region_file (file) : Region file
        output_file (string) : Output file name
        macs_summits_file (file) : MACS summits file
    """
    #initialize chromosomal positions based on gff provided
    opt_chrom, opt_start, opt_stop, opt_peakid, opt_score, opt_none, opt_gene, opt_transcript = 0, 1, 2, 3, 4, 6, 8, 9
    summits_dict = {}
    if macs_summits_file: #read summits file
        with open(macs_summits_file, 'r') as summits_input:
            summits_dict = {line.strip().split('\t')[opt_peakid]: (line.strip().split('\t')[opt_start], line.strip().split('\t')[opt_stop]) for line in summits_input}

    #read intersectbed file
    with open(region_file, 'r') as intersect_input, open(output_file, 'w') as final_output:
        peaks_dict, peaks_score, gene_dict, transcript_dict, peak_region = {}, {}, {}, {}, []
        for line in intersect_input:
            line = line.strip().split('\t')
            if int(line[opt_none]) != -1:
                string_identifier = tuple(line[i] for i in [opt_chrom, opt_start, opt_stop])
                peaks_dict[string_identifier] = line[opt_peakid]
                peaks_score[line[opt_peakid]] = line[opt_score]
                peak_region.append(string_identifier)
                gene_dict[string_identifier] = gene_dict.get(string_identifier, "") + line[opt_gene]
                transcript_dict[string_identifier] = transcript_dict.get(string_identifier, "") + line[opt_transcript]

        unique_peak_region = list(set(peak_region))

        header = "peak_id\tchrom\tpeak_start\tpeak_end\tsummit_start\tsummit_end\tgene_name\ttranscript_id\tpeak_score\n" if macs_summits_file else "peak_id\tchrom\tpeak_start\tpeak_end\tgene_name\ttranscript_id\tpeak_score\n"

        final_output.write(header)
        for chromregion in unique_peak_region:
            peakid = peaks_dict[chromregion]
            if macs_summits_file:
                final_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(peakid, '\t'.join(chromregion), '\t'.join(summits_dict[peakid]),','.join(unique(gene_dict[chromregion])),','.join(unique(transcript_dict[chromregion])),peaks_score[peakid]))
            else:
                final_output.write("%s\t%s\t%s\t%s\t%s\n" %(peakid, '\t'.join(chromregion), ','.join(unique(gene_dict[chromregion])), ','.join(unique(transcript_dict[chromregion])), peaks_score[peakid]))


def unique(list1):
    """ Get unique list
    Args:
        list1 (list) : list of variable separated by comma ','

    """
    return list(set(list1.split(",")))


def include_genes(macs_peaks_file, file_name):
    """ Annotate Genes in peak regions
    Args:
        macs_peaks_file (file) : MACS peaks file
        output_file (string) : Output file name
    """
    peaks_chr_dict, peaks_count, gene_names_dict, region_dict, peaks_id_dict, region_names_list, master_generegion_dict, barcoderegions_dict, barcodepeaks_id_dict, master_peaksid_dict = {}, 0, {}, defaultdict(dict), defaultdict(dict), [], defaultdict(dict), {}, {}, defaultdict(dict)

    #reading peaks
    peaks_chr_dict = {line.split('\t')[3]: line.split('\t')[0:3] for line in open(macs_peaks_file, 'r')}
    peaks_count = sum(1 for _ in open(macs_peaks_file, 'r'))

    #reading gtf gene names
    gene_names_dict = {(line.split('\t')[0],line.split('\t')[1]): line.split('\t')[0] for line in open("tempPrefix.genes.names", 'r')}

    #reading the genomic regions files
    for region_prefix in ("promoter", "genebody", "window", "closest"):
        region = region_prefix + ".tempPrefix"
        region_name = Path(region).stem.upper()
        region_names_list.append(region_name)
        with open(region, 'r') as region_input:
            for line in region_input:
                line = line.strip().split('\t')
                region_dict[region_name][(line[8], line[9])] = line[8:10]
                peaks_id_dict[region_name][line[3]] = line[0:3]
    peaks_id_list = peaks_chr_dict.keys()
    peaks_id_list = sorted(peaks_id_list, key=lambda item:
                           (int(item[10:]) if item[10:].isdigit()
                            else print(item[10:], type(item[10:])), item))

    #order regions and gene names
    region_names_list.sort()
    orderedgene_names_dict = collections.OrderedDict(sorted(gene_names_dict.items()))

    #barcoding the identified regions
    for eachgene in orderedgene_names_dict:
        for eachkey in region_names_list:
            found = 1 if eachgene in region_dict[eachkey] else 0
            master_generegion_dict[eachkey][eachgene] = found
            barcoderegions_dict[eachgene] = barcoderegions_dict.get(eachgene, "") + str(found)

    for eachpeaks_id in peaks_id_list:
        for eachkey in region_names_list:
            found = 1 if eachpeaks_id in peaks_id_dict[eachkey] else 0
            master_peaksid_dict[eachkey][eachpeaks_id] = found
            barcodepeaks_id_dict[eachpeaks_id] = barcodepeaks_id_dict.get(eachpeaks_id, "") + str(found)

    
    #write to output file
    gene_file_name = file_name.replace('.txt', '.genes.txt')

    with open(gene_file_name, 'w') as final_output:
        final_output.write("GeneName\tTranscript_ID\tBARCODE\t%s\tEVER\n" % ('\t'.join(region_names_list)))
        for eachgene in orderedgene_names_dict:
            string_identifier = [barcoderegions_dict[eachgene]]
            found_sum = sum(master_generegion_dict[region][eachgene] for region in region_names_list)
            string_identifier.extend(str(master_generegion_dict[region][eachgene]) for region in region_names_list)
            string_identifier.append("1" if found_sum >= 1 else "0")

            final_output.write("%s\t%s\n" % ('\t'.join(eachgene), '\t'.join(string_identifier)))

    peaks_file_name = file_name.replace('.txt', '.peaks.txt')
    pdf_file_name = file_name.replace('.txt', '.distribution.pdf')
    with open(peaks_file_name, 'w') as final_output_peaks:
        final_output_peaks.write("chrom\tstart\tstop\tpeak_id\tBARCODE\t%s\n" %('\t'.join(region_names_list)))
        for eachpeaks_id in peaks_id_list:
            string_identifier = [barcodepeaks_id_dict[eachpeaks_id]]
            string_identifier.extend(str(master_peaksid_dict[region][eachpeaks_id]) for region in region_names_list)
            final_output_peaks.write("%s\t%s\t%s\n" %('\t'.join(peaks_chr_dict[eachpeaks_id]), eachpeaks_id, '\t'.join(string_identifier)))

    region_count_dict = {region: round(len(peaks_id_dict[region]) * 100 / peaks_count) for region in ("PROMOTER", "GENEBODY", "WINDOW")}
    rects = plt.bar(region_count_dict.keys(), region_count_dict.values(), width=0.4)
    axes = plt.gca()
    axes.set_ylim([0, 100])
    peaks_title = Path(macs_peaks_file).stem
    plt.title(peaks_title.split('.sorted')[0])
    plt.ylabel("% of peaks")

    #label the bars
    xpos = 'center'
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                 '{}%'.format(height), ha=ha[xpos], va='bottom')

    #save figure
    plt.savefig(pdf_file_name)

    for eachregion in ("PROMOTER", "GENEBODY", "WINDOW"):
        print("Total number of %s = %d" %(eachregion, len(peaks_id_dict[eachregion])))
        print("Percentage of %s = %.3f" %(eachregion, len(peaks_id_dict[eachregion])*100/peaks_count))


def main():
    '''
    Generate genomic coordinates of all promoters, window, genebody, TSS

    Definitions: promoters: +/- 1kb of TSS (transcription start site)
                 window: 10kb upstream --> 3kb downstream of genebody
                 genebody: 1kb upstream of TSS --> TTS/TES (transcription termination/end site)
                 TSS: TSS (transcription start site)

    Annotate Peaks based on genomic features listed above

    Comparison master files for all genes/transcripts
        and their genomic regions within peaks called.
    Genomic regions are specified based on regions provided,
        default expectation are promoter, genebody, window, closest

    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", dest="gtf", required=True, help="Enter .gtf/gff file to be processed.")
    parser.add_argument("-c", "--chrom", dest="chrom_sizes", required=True, help="Enter ucsc chrom sizes file to be processed.")
    parser.add_argument("-p", "--peaks", dest="macs_peaks", required=True, help="Enter MACS peaks bed file.")
    parser.add_argument("-s", "--summit", dest="macs_summits", required=False, help="Enter MACS summit bed file.")

    options = parser.parse_args()

    # Stage 1
    # Extract genomic regions of interest
    gtf_to_genes(options.chrom_sizes, options.gtf)

    #Peaks filecheck
    checkcolumns = sum(1 for _ in open(options.macs_peaks, 'r'))
    if checkcolumns <= 0: sys.exit("ERROR :\tNo peaks were identified\n") 

    # Turn 4 column bed to 5 (for SICER peaks: include unique id)
    with open(options.macs_peaks, 'r') as peaks_file:
        numberofcolumns = len(peaks_file.readline().strip().split('\t'))
        fifthcol_value = next(peaks_file).split('\t')[4].strip() #FDRisland from SICER produce a fifth empty column
        
    new_macs_peaks = options.macs_peaks

    if numberofcolumns == 4 or len(fifthcol_value) < 1:
        new_macs_peaks = options.macs_peaks + ".tempPrefix"
        # command = "cat " + options.macs_peaks + ' | awk -F\\\\t ' + "'" + \
        #           '{print $1 "\\t" $2 "\\t" $3 "\\t" $1":"$2"-"$3 "\\t" $4}' + \
        #           "' > tempPrefix-new.macs_peaks.bed"
        # print(command)
        # os.system(command)
        macs_peaks_output = open(new_macs_peaks, 'w')
        macs_peaks_input = open(options.macs_peaks,'r')
        macs_peaks_count = 1
        for line in macs_peaks_input:
            line = line.strip().split('\t')
            new_peaks_id = "SICERpeak_" + str(macs_peaks_count)
            results = ("{0}\t{1}\t{2}\n".format("\t".join(line[0:3]), new_peaks_id, line[3]))
            macs_peaks_output.write(results)
            macs_peaks_count += 1
        macs_peaks_output.close()
        
    # Convert TSS to compatible bedtools (closestBed) format and obtain center of peaks coordinate
    command = "cut -f1,4,5,9,10 annotation/TSS-region_genes.gff | sort -k1,1 -k2,2n " + \
              "> tempPrefix-sorted.TSS; " + \
              "cat " + new_macs_peaks + ' | awk -F\\\\t ' + "'" + \
              '{print $1 "\\t" int(($2+$3)/2) "\\t" int(($2+$3)/2) "\\t" $4 "\\t" $5}' + \
              "' | sort -k1,1 -k2,2n > tempPrefix-centerofpeaks.bed;" + \
              "closestBed -t all -a tempPrefix-centerofpeaks.bed -b tempPrefix-sorted.TSS" + \
              "> closest.tempPrefix"
    os.system(command)

    # Stage 2
    # Annotate peaks based on previously extracted genomic coordinates
    for region in ("promoter", "genebody", "window"):
        command = 'intersectBed -b ' + new_macs_peaks + ' -a annotation/' + region + \
                  '-region_genes.gff -wb | cut -f1,4,5,9,10,11,12,13,14,15 |' + "awk -F\\\\t '" + \
                  '{print $6 "\\t" $7 "\\t" $8 "\\t" $9 "\\t" $10 "\\t" $1 "\\t" $2 "\\t" $3 ' + \
                  '"\\t" $4 "\\t" $5}' + \
                  "' > " + region + ".tempPrefix"
        os.system(command)

    # Convert annotated regions to tab delimited file and include summit location
    for region in ("promoter", "genebody", "window", "closest"):
        region_file = region + ".tempPrefix"
        output_file = "peaks_within_" + region + ".regions.txt"
        if options.macs_summits:
            annotate_regions(region_file, output_file, options.macs_summits)
        else:
            annotate_regions(region_file, output_file)

    command = "mv peaks_within_closest.regions.txt centerofpeaks_closest.regions.txt; \
              cut -f 9,10 annotation/original_genes.gff > tempPrefix.genes.names"
    os.system(command)

    # Stage 3
    # Binary matrix of all genes
    include_genes(new_macs_peaks, "peaks_compared_regions.txt")

    command = "rm -rf *tempPrefix*"
    os.system(command)


if __name__ == "__main__":
    main()
