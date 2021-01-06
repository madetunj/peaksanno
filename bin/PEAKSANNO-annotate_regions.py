#!/usr/bin/env python3

'''
    Annotate Peaks based on genomic features provided
'''

import argparse

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
    parser.add_argument("-g", "--gff", dest="intersect_file", required=True,
                        help="Enter 'intersectBed -wb' file generated.")
    parser.add_argument("-s", "--summit", dest="macs_summits_file", required=True,
                        help="Enter MACS summit bed file.")
    parser.add_argument("-o", "--output", dest="output_file",
                        default="peaks_within.regions.txt", help="Enter output file name")

    options = parser.parse_args()

    #initialize chromosomal positions based on gff provided
    opt_chrom = 0
    opt_start = 1
    opt_stop = 2
    opt_peakid = 3
    opt_score = 4
    opt_gene = 8
    opt_transcript = 9

    #read summits file
    summits_input = open(options.macs_summits_file, 'r')
    summits_dict = {}
    for line in summits_input:
        line = line.strip('\n').split('\t')
        string_identifier = line[opt_start], line[opt_stop]
        summits_dict[line[opt_peakid]] = string_identifier

    #read intersectbed file
    intersect_input = open(options.intersect_file, 'r')
    peaks_dict = {}
    peaks_score = {}
    gene_dict = {}
    transcript_dict = {}
    peak_region = []
    for line in intersect_input:
        line = line.strip('\n').split('\t')
        string_identifier = line[opt_chrom], line[opt_start], line[opt_stop]
        peaks_dict[string_identifier] = line[opt_peakid]
        peaks_score[line[opt_peakid]] = line[opt_score]
        peak_region.append(string_identifier)

        if string_identifier in gene_dict:
            gene_dict[string_identifier] = ("{0},{1}".format(gene_dict[string_identifier],
                                                             line[opt_gene]))
            transcript_dict[string_identifier] = ("{0},{1}".
                                                  format(transcript_dict[string_identifier],
                                                         line[opt_transcript]))
        else:
            gene_dict[string_identifier] = line[opt_gene]
            transcript_dict[string_identifier] = line[opt_transcript]

    #each peak id is unique
    unique_peak_region = []
    for  eachpeak_region in peak_region:
        if eachpeak_region not in unique_peak_region:
            unique_peak_region.append(eachpeak_region)

    #write to output file
    final_output = open(options.output_file, 'w')
    final_output.write("peak_id\tchrom\tpeak_start\tpeak_end\t"
                       "summit_start\tsummit_end\tgene_name\t"
                       "transcript_id\tpeak_score\n")

    #print regions identified
    for chromregion in unique_peak_region:
        peakid = peaks_dict[chromregion]
        final_output.write("%s\t%s\t%s\t%s\t%s\t%s\n"
                           %(peakid, '\t'.join(chromregion), '\t'.join(summits_dict[peakid]),
                             ','.join(unique(gene_dict[chromregion])),
                             ','.join(unique(transcript_dict[chromregion])),
                             peaks_score[peakid]))

if __name__ == "__main__":
    main()
