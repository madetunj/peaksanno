#!/usr/bin/env python3
'''
    Comparison master files for all genes/transcripts
        and their genomic regions within peaks called.
    Genomic regions are specified based on regions provided,
        default expectation are promoter, genebody, window, closest

'''

import argparse
from pathlib import Path
import collections
from collections import defaultdict

def unique(list1):
    """ Get unique list
    Args:
        list1 (list) : list of variable separated by pipes'||'

    """

    # intilize a null list
    unique_list = []

    # traverse for all elements
    for variable in list1.split("||"):
        if variable not in unique_list:
            unique_list.append(variable)
    return unique_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genes", dest="gene_names_file", required=True,
                        help="Enter gene names file")
    parser.add_argument("-f", "--files", dest="regions_files", required=True,
                        nargs='+', help="Enter 1 or more genomic regions file(s).")
    parser.add_argument("-o", "--output", dest="output_file",
                        default="peaks_comparison.txt", help="Enter output file name")

    options = parser.parse_args()

    #reading the main genes file
    gene_names_dict = {}
    gene_names_input = open(options.gene_names_file, 'r')
    for line in gene_names_input:
        line = line.strip('\n').split('\t')
        gene_names_dict[(line[0], line[1])] = line[0]

    #reading the genomic regions files
    region_names_list = []
    region_dict = defaultdict(dict)
    for region in options.regions_files:
        region_name = Path(region).stem
        region_names_list.append(region_name.upper())
        region_input = open(region, 'r')
        for line in region_input:
            line = line.strip('\n').split('\t')
            region_dict[region_name.upper()][(line[8], line[9])] = (line[8], line[9])

    #order regions and gene names
    region_names_list.sort()
    orderedgene_names_dict = collections.OrderedDict(sorted(gene_names_dict.items()))

    barcoderegions_dict = {}
    master_generegion_dict = defaultdict(dict)
    for eachgene in orderedgene_names_dict:
        for eachkey in region_names_list:
            found = 0
            if eachgene in region_dict[eachkey]:
                found = 1
            master_generegion_dict[eachkey][eachgene] = found
            if eachgene in barcoderegions_dict:
                barcoderegions_dict[eachgene] += str(found)
            else:
                barcoderegions_dict[eachgene] = str(found)

    #write to output file
    final_output = open(options.output_file, 'w')
    final_output.write("GeneName\tTranscript_ID\tBARCODE\t%s\tEVER\n"
                       %('\t'.join(region_names_list)))

    for eachgene in orderedgene_names_dict:
        string_identifier = [barcoderegions_dict[eachgene]]
        found_sum = 0
        for eachkey in region_names_list:
            string_identifier.append(str(master_generegion_dict[eachkey][eachgene]))
            found_sum += master_generegion_dict[eachkey][eachgene]
        if found_sum >= 1:
            string_identifier.append("1")
        else:
            string_identifier.append("0")

        final_output.write("%s\t%s\n" %('\t'.join(eachgene), '\t'.join(string_identifier)))

if __name__ == "__main__":
    main()
