#!/usr/bin/env python3

'''
    Generate genomic coordinates of all promoters, window, genebody, TSS

    Definitions: promoters: +/- 1kb of TSS (transcription start site)
                 window: 10kb upstream --> 3kb downstream of genebody
                 genebody: 1kb upstream of TSS --> TTS/TES (transcription termination/end site)
                 TSS: TSS (transcription start site)

'''

import os
import sys
import re
import argparse

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", dest="gtf_file", required=True,
                        help="Enter .gtf/gff file to be processed.")
    parser.add_argument("-c", "--chrom", dest="chrom_sizes", required=True,
                        help="Enter ucsc chrom sizes file to be processed.")

    options = parser.parse_args()
    #print(options)

    #initialize outputfiles
    original_output = open('annotation/original_genes.gff', 'w')
    promoters_output = open('annotation/promoter-region_genes.gff', 'w')
    window_output = open('annotation/window-region_genes.gff', 'w')
    tss_output = open('annotation/TSS-region_genes.gff', 'w')
    genebody_output = open('annotation/genebody-region_genes.gff', 'w')

    #read ucsc chrom sizes, gtf
    chrom_sizes_input = open(options.chrom_sizes, 'r')
    chrom_sizes_dict = {}
    for line in chrom_sizes_input:
        lines = line.split('\t')
        chrom_sizes_dict[lines[0]] = lines[1].rstrip("\n")

    #feature = options.feature
    feature_dict = {}
    if not options.gtf_file.split('.')[-1] == 'gtf':
        sys.exit("ERROR :\tGTF required")

    #reading the feature type to be used
    gtf_input = open(options.gtf_file, 'r')
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
    gtf_input = open(options.gtf_file, 'r')
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


if __name__ == "__main__":
    main()
