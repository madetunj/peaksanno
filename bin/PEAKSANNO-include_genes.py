#!/usr/bin/env python3
'''
    Comparison master files for all genes/transcripts
        and their genomic regions within peaks called.
    Genomic regions are specified based on regions provided,
        default expectation are promoter, genebody, window, closest

'''

import argparse
import matplotlib.pyplot as plt
from pathlib import Path
import collections
from collections import defaultdict
     
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genes", dest="gene_names_file", required=True,
                        help="Enter gene names file")
    parser.add_argument("-p", "--peaks", dest="peaks_bed_file", required=True,
                        help="Enter MACS peaks bed file")
    parser.add_argument("-f", "--files", dest="regions_files", required=True,
                        nargs='+', help="Enter 1 or more genomic regions file(s).")

    options = parser.parse_args()

    #counting total number of peaks
    file_name = "peaks_compared_regions.txt"
    
    peaks_input = open(options.peaks_bed_file,'r')
    peaks_chr_dict = {}
    peaks_count = 0
    for line in peaks_input:
        line = line.strip('\n').split('\t')
        peaks_count += 1
        peaks_chr_dict[line[3]] = (line[0:3])
        

    #reading the gtf genes file
    gene_names_dict = {}
    gene_names_input = open(options.gene_names_file, 'r')
    for line in gene_names_input:
        line = line.strip('\n').split('\t')
        gene_names_dict[(line[0], line[1])] = line[0]

    #reading the genomic regions files
    region_names_list = []
    region_dict = defaultdict(dict)
    peaks_id_dict = defaultdict(dict)
    for region in options.regions_files:
        region_name = Path(region).stem
        region_names_list.append(region_name.upper())
        region_input = open(region, 'r')
        for line in region_input:
            line = line.strip('\n').split('\t')
            region_dict[region_name.upper()][(line[8], line[9])] = (line[8], line[9])
            peaks_id_dict[region_name.upper()][line[3]] = (line[0:3])
    peaks_id_list = peaks_chr_dict.keys()
    peaks_id_list = sorted(peaks_id_list, key=lambda item: (int(item[10:]) if item[10:].isdigit() else print(item[10:],type(item[10:])), item))

    #order regions and gene names
    region_names_list.sort()
    orderedgene_names_dict = collections.OrderedDict(sorted(gene_names_dict.items()))

    #barcoding the identified regions
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

    barcodepeaks_id_dict = {}
    master_peaksid_dict = defaultdict(dict)
    for eachpeaks_id in peaks_id_list:
        for eachkey in region_names_list:
            found = 0
            #print(eachpeaks_id)
            if eachpeaks_id in peaks_id_dict[eachkey]:
                found = 1
                #print("yes")
            master_peaksid_dict[eachkey][eachpeaks_id] = found
            if eachpeaks_id in barcodepeaks_id_dict:
                barcodepeaks_id_dict[eachpeaks_id] += str(found)
            else:
                barcodepeaks_id_dict[eachpeaks_id] = str(found)
   
    #write to output file
    gene_file_name = file_name.replace('.txt','.genes.txt')
    
    final_output = open(gene_file_name, 'w')
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
        
    peaks_file_name = file_name.replace('.txt','.peaks.txt')
    pdf_file_name = file_name.replace('.txt','.distribution.pdf')
    final_output_peaks = open(peaks_file_name, 'w')
    final_output_peaks.write("chrom\tstart\tstop\tpeak_id\tBARCODE\t%s\n"
                       %('\t'.join(region_names_list)))
    for eachpeaks_id in peaks_id_list:
        string_identifier = [barcodepeaks_id_dict[eachpeaks_id]]
        for eachkey in region_names_list:
            string_identifier.append(str(master_peaksid_dict[eachkey][eachpeaks_id]))

        final_output_peaks.write("%s\t%s\t%s\n" %('\t'.join(peaks_chr_dict[eachpeaks_id]), eachpeaks_id, '\t'.join(string_identifier)))
    
    
    #the bar plot
    region_count_dict = {}    
    #for eachregion in region_count_dict:
    region_count_dict["Promoter\n(1kb up/down TSS)"] = round(len(peaks_id_dict["PROMOTER"])*100/peaks_count)
    region_count_dict["Genebody\n(1kb up TSS/TES)"] = round(len(peaks_id_dict["GENEBODY"])*100/peaks_count)
    region_count_dict["Window\n(10kb up TSS/3kb down TES)"] = round(len(peaks_id_dict["WINDOW"])*100/peaks_count)
    
    rects = plt.bar(region_count_dict.keys(),region_count_dict.values(), width=0.4)
    axes = plt.gca()
    axes.set_ylim([0,100])
    peaks_title=Path(options.peaks_bed_file).stem
    plt.title(peaks_title.split('.sorted')[0])
    plt.ylabel("% of peaks")

    #label the bars
    xpos='center'
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off

    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,
                '{}%'.format(height), ha=ha[xpos], va='bottom')
        
    #save figure
    plt.savefig(pdf_file_name)

    for eachregion in ("PROMOTER","GENEBODY","WINDOW"):
        print("Total number of %s = %d" %(eachregion,len(peaks_id_dict[eachregion])))
        print ("Percentage of %s = %.3f" %(eachregion,len(peaks_id_dict[eachregion])*100/peaks_count))

if __name__ == "__main__":
    main()
