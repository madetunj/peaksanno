#!/usr/bin/env python3

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
    
import os
import sys
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
    lines[3] = int(lines[3])
    lines[4] = int(lines[4])

    if lines[6] == "+":
        if tss:
            end_region = 4
        else:
            end_region = 3

        end_position = lines[end_region] + down_dist
        start_position = lines[3] - up_dist

    elif lines[6] == "-":
        if tss:
            end_region = 3
        else:
            end_region = 4

        start_position = lines[end_region] - down_dist
        end_position = lines[4] + up_dist

    if start_position < 1:
        start_position = 1

    if end_position > int(chromz[lines[0]]):
        end_position = chromz[lines[0]]

    output_gff = ("{0}\t{1}\t{2}\t{3}\n".format("\t".join(lines[0:3]),
                                                start_position, end_position, "\t".join(lines[5:])))
    return output_gff


def gtf_to_genes(chrom_sizes_file, gtf_file):
    """ GTF to GENES
    Args:
        chrom_sizes_file (file) : Chromosomal sizes file
        gtf_file (file) : GTF file
    """
    #initialize outputfiles
    original_output = open('annotation/original_genes.gff', 'w')
    promoters_output = open('annotation/promoter-region_genes.gff', 'w')
    window_output = open('annotation/window-region_genes.gff', 'w')
    tss_output = open('annotation/TSS-region_genes.gff', 'w')
    genebody_output = open('annotation/genebody-region_genes.gff', 'w')

    #read ucsc chrom sizes, gtf
    chrom_sizes_input = open(chrom_sizes_file, 'r')
    chrom_sizes_dict = {}
    for line in chrom_sizes_input:
        lines = line.split('\t')
        chrom_sizes_dict[lines[0]] = lines[1].rstrip("\n")

    #feature = options.feature
    feature_dict = {}
    if not gtf_file.split('.')[-1] == 'gtf':
        sys.exit("ERROR :\tGTF required")

    #reading the feature type to be used
    gtf_input = open(gtf_file, 'r')
    for line in gtf_input:
        if not line.startswith('#'):
            lines = line.split("\t")
            feature_dict[lines[2]] = lines[2]

    if 'transcript' in feature_dict:
        feature = "transcript"
    elif 'gene' in feature_dict:
        feature = "gene"
    else:
        sys.exit("ERROR :\tGTF with either transcript/gene annotation is needed")
    print("NOTE :\tFeature type used is '%s'" %(feature))

    #organize gtf file
    gtf_input = open(gtf_file, 'r')
    for line in gtf_input:
        if not line.startswith('#'):
            lines = line.split("\t")
            newline = lines[8].split(' ')
            gene_name = re.sub('[\"\;]', '', newline[1]) #clean gene_name
            transcript_id = re.sub('[\"\;]', '', newline[3]) #clean transcript_id
            if lines[2] == feature:
                if re.match("chr", lines[0]):
                    new_chr_name = lines[0]
                else:
                    new_chr_name = "chr"+lines[0]

                results = ("{0}\t{1}\t{2}\t{3}".format(new_chr_name,
                                                       "\t".join(lines[1:8]),
                                                       gene_name, transcript_id))

                #extract chromosomal coordinates and store in outputfiles
                original_output.write(results + "\n")
                final_output = parse_genelocations(chrom_sizes_dict, results, 1000, 1000, False)
                promoters_output.write(final_output)
                final_output = parse_genelocations(chrom_sizes_dict, results, 10000, 3000, True)
                window_output.write(final_output)
                final_output = parse_genelocations(chrom_sizes_dict, results, 0, 0, False)
                tss_output.write(final_output)
                final_output = parse_genelocations(chrom_sizes_dict, results, 1000, 0, True)
                genebody_output.write(final_output)


def unique(list1):
    """ Get unique list
    Args:
        list1 (list) : list of variable separated by comma ','

    """

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for variable in list1.split(","):
        if variable not in unique_list:
            unique_list.append(variable)
    return unique_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", dest="gtf", required=True,
                        help="Enter .gtf/gff file to be processed.")
    parser.add_argument("-c", "--chrom", dest="chrom_sizes", required=True,
                        help="Enter ucsc chrom sizes file to be processed.")
    
    #parser.add_argument("-n", "--intersect", dest="intersect", required=True,
    #                    help="Enter 'intersectBed -wb' file generated.")
    parser.add_argument("-s", "--summit", dest="macs_summits", required=True,
                        help="Enter MACS summit bed file.")
    parser.add_argument("-o", "--output", dest="output",
                        default="peaks_within.regions.txt", help="Enter output file name")
    
    # parser.add_argument("-g", "--genes", dest="gene_names_file", required=True,
    #                     help="Enter gene names file")
    # parser.add_argument("-p", "--peaks", dest="peaks_bed_file", required=True,
    #                     help="Enter MACS peaks bed file")
    # parser.add_argument("-f", "--files", dest="regions_files", required=True,
    #                     nargs='+', help="Enter 1 or more genomic regions file(s).")
    # 

    options = parser.parse_args()
    #print(options)
    
    gtf_to_genes(options.chrom_sizes, options.gtf)
    
    for region in ("promoter","genebody","window"):
        command = 'intersectBed -b $BEDFILE -a annotation/' + region + '-region_genes.gff -wb | cut -f1,4,5,9,10,11,12,13,14,15 |' + "awk -F\\\\t '" + '{print $6 "\\t" $7 "\\t" $8 "\\t" $9 "\\t" $10 "\\t" $1 "\\t" $2 "\\t" $3 "\\t" $4 "\\t" $5}' + "' > " + region + ".tempPrefix"
        print (command)
        stats = subprocess.Popen(command,stdin = subprocess.PIPE,stderr = subprocess.PIPE,stdout = subprocess.PIPE,shell = True)
        statLines = stats.stdout.readlines()
        stats.stdout.close()
        print (statLines)
#     command = 'intersectBed -b $BEDFILE -a annotation/%s-region_genes.gff -wb
#     annotationArray=( "promoter" "genebody" "window" )
# for args in "${annotationArray[@]}"
# do
#   # intersect similar peaks with gene regions
#   intersectBed -b $BEDFILE -a annotation/$args*-region_genes.gff -wb | \
#   cut -f1,4,5,9,10,11,12,13,14,15 | \
#   awk -F\\t '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' \
#   > $args.$tempPrefix
#     #initialize chromosomal positions based on gff provided
#     opt_chrom = 0
#     opt_start = 1
#     opt_stop = 2
#     opt_peakid = 3
#     opt_score = 4
#     opt_gene = 8
#     opt_transcript = 9
# 
#     #read summits file
#     summits_input = open(options.macs_summits, 'r')
#     summits_dict = {}
#     for line in summits_input:
#         line = line.strip('\n').split('\t')
#         string_identifier = line[opt_start], line[opt_stop]
#         summits_dict[line[opt_peakid]] = string_identifier
# 
#     #read intersectbed file
#     intersect_input = open(options.intersect, 'r')
#     peaks_dict = {}
#     peaks_score = {}
#     gene_dict = {}
#     transcript_dict = {}
#     peak_region = []
#     for line in intersect_input:
#         line = line.strip('\n').split('\t')
#         string_identifier = line[opt_chrom], line[opt_start], line[opt_stop]
#         peaks_dict[string_identifier] = line[opt_peakid]
#         peaks_score[line[opt_peakid]] = line[opt_score]
#         peak_region.append(string_identifier)
# 
#         if string_identifier in gene_dict:
#             gene_dict[string_identifier] = ("{0},{1}".format(gene_dict[string_identifier],
#                                                              line[opt_gene]))
#             transcript_dict[string_identifier] = ("{0},{1}".
#                                                   format(transcript_dict[string_identifier],
#                                                          line[opt_transcript]))
#         else:
#             gene_dict[string_identifier] = line[opt_gene]
#             transcript_dict[string_identifier] = line[opt_transcript]
# 
#     #each peak id is unique
#     unique_peak_region = []
#     for  eachpeak_region in peak_region:
#         if eachpeak_region not in unique_peak_region:
#             unique_peak_region.append(eachpeak_region)
# 
#     #write to output file
#     final_output = open(options.output, 'w')
#     final_output.write("peak_id\tchrom\tpeak_start\tpeak_end\t"
#                        "summit_start\tsummit_end\tgene_name\t"
#                        "transcript_id\tpeak_score\n")
# 
#     #print regions identified
#     for chromregion in unique_peak_region:
#         peakid = peaks_dict[chromregion]
#         final_output.write("%s\t%s\t%s\t%s\t%s\t%s\n"
#                            %(peakid, '\t'.join(chromregion), '\t'.join(summits_dict[peakid]),
#                              ','.join(unique(gene_dict[chromregion])),
#                              ','.join(unique(transcript_dict[chromregion])),
#                              peaks_score[peakid]))
# 


if __name__ == "__main__":
    main()
