#!/usr/bin/env python
# coding=utf-8

"""
This script is to load the config file for nanostr, including:
1. locus information: locus, position, prefix seq, suffix seq
2. struction info: fact-sheet
3. reads from bam file using samtools with flag 2308(unmapped, not pri, sup)
"""

import re, os
from random import seed, sample
from collections import defaultdict

## locus information: locus, position, prefix seq, suffix seq
def locus_info_loader(locus_info_file):
    locus_info_dict = defaultdict(list)

    with open(locus_info_file, 'r') as file:
        for ind, line in enumerate(file):
            if ind == 0:
                continue
            if line.strip() == "":
                continue
            ln_lst   = line.strip().split(",")
            locus, chrom, start, end, flank_len, prefix_seq, suffix_seq = \
            ln_lst[0], ln_lst[1], ln_lst[2], ln_lst[3], ln_lst[4], ln_lst[5], ln_lst[6]
            enlarged = "%s:%s-%s" % (chrom, int(start) - int(flank_len), int(end) + int(flank_len))
            locus_info_dict[locus] = [enlarged, flank_len, prefix_seq, suffix_seq]

    return locus_info_dict


## struction info: fact-sheet
def pattern_loader(fact_sheet_path, locus):
    def fact_line_finder(line):
        counted_pat_dict = {}
        full_pat_dict    = {}

        for i in re.findall('\[([A-Z]+)\]', line):
            if i.upper() not in counted_pat_dict:
                counted_pat_dict[i] = 'X%s' % len(counted_pat_dict)
                full_pat_dict[i]    = 'X%s' % len(full_pat_dict)

        for i in re.findall('\s([a-z]+)', line):
            if i.upper() not in full_pat_dict:
                full_pat_dict[i] = 'x%s' % len(full_pat_dict)

        ordered_pat_list = re.findall('\s+\[*([ATCGatcg]+)', line)

        return counted_pat_dict, full_pat_dict, ordered_pat_list
    
    counted_pat_dict, full_pat_dict, ordered_pat_list = None, None, None
    with open(fact_sheet_path, 'r') as file:
        for line in file:
            if re.findall(locus.upper(), line.upper()):
                counted_pat_dict, full_pat_dict, ordered_pat_list = fact_line_finder(line)
    
    if counted_pat_dict is None:
        return None, None, None
    
    ordered_pat_deno_list = None
    # ordered_pat_deno_list = [ full_pat_dict[i] for i in ordered_pat_list ]
    # not_counted_unit_list = [ i for i in ordered_pat_list if i not in counted_pat_dict]
    # not_counted_deno_list = [ full_pat_dict[i.upper()] for i in ordered_pat_list if i not in counted_pat_dict]

    return counted_pat_dict, ordered_pat_deno_list, ordered_pat_list # , not_counted_unit_list, not_counted_deno_list


## reads
def bam_loader(samtools, bamfile_path, enlarge_region, depth):
    bamin   = os.popen('%s view -F 2308 %s %s' % (samtools, bamfile_path, enlarge_region))
    bam_lst = bamin.readlines()
    
    readid = []
    reads  = []

    for cnt, line in enumerate(bam_lst):
        lineLst = line.strip().split()
        readid.append(lineLst[0])
        reads.append(lineLst[9])

    if len(reads) <= depth or depth == -1:
        return readid, reads
    else:
        seed(1)
        index_lst = list(range(len(readid)))
        choosen_lst = sample(index_lst, depth)
        downsampled_readid = [readid[i] for i in choosen_lst]
        downsampled_read   = [reads[i] for i in choosen_lst]
        return downsampled_readid, downsampled_read