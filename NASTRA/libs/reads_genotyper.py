#!/usr/bin/env python
# coding=utf-8
from re import findall
# from .raw_reads_processor import aln_score_func


def reads_genotyper(locus, bracket_allele):
    
    unit_cnt_lst = [float(i) for i in findall('\d+[.\d]*', bracket_allele)]
    if unit_cnt_lst == []:
        genotype = 0
    else:
        genotype = sum(unit_cnt_lst)


    # if locus in ['D20S482']:
    #     return bracket_allele, genotype

    ## manully editted

    if locus.upper() == 'DXS10135':  # calculated by sequence length
        seq_len = 0
        for unit in bracket_allele.split():
            if '[' not in unit:
                seq_len += len(unit.strip())
            else:
                motif    = findall('\[([A-Za-z]*)\]', unit)[0]
                count    = findall('\d+', unit)[0]
                seq_len += len(motif) * int(count)

        genotype = int( (seq_len - 13) / 4 ) #  -7 GAAAGGA; ignore head and tail 
        return bracket_allele, genotype
    
    if locus.upper() == 'DXS10103':
        seq_lst       = bracket_allele.split(' ')[1:-1]
        div = sum([len(i) for i in seq_lst if '[' not in i]) // 4
        mod = sum([len(i) for i in seq_lst if '[' not in i]) % 4
        genotype += div + mod / 10
       
        genotype -= 1
        return bracket_allele, genotype

    if locus.upper() == "DXS10074":
        cnts = findall("\](\d+)", bracket_allele)
        ## actually, genotype: AAA [AAGA]1 AGA [AAGA]13 [AAGG]1 [AAGA]
        ## equals to AAAAAG [AAGA]14 [AAGG]1 [AAGA]1
        genotype = sum([int(i) for i in cnts])
        return bracket_allele, genotype

    
    ## genotype automatically
    # editted_allele = allele_corrector(bracket_allele, locus)
    editted_allele = bracket_allele 
    genotype       = single_allele_genotyper(locus, editted_allele)
    return editted_allele, genotype




def single_allele_genotyper(locus, editted_allele):

    if editted_allele == 'others':
        return 'others'

    unit_cnt_lst = [float(i) for i in findall('\d+[.\d]*', editted_allele)]

    if unit_cnt_lst == []:
        allele = 0
    else:
        allele = sum(unit_cnt_lst)

    ## TODO: fixed deviation to revise
    if locus == "D13S317":
        seq_lst       = editted_allele.split(' ')
        deviation0    = abs( len(seq_lst[-1].strip()) - len("AATCAATCATC") - 3 )      # tail with 3 more nucs
        deviation1    = abs( len(seq_lst[-1].strip()) - len("AATCATC") - 3 )          # tail with 3 more nucs
        deviation2    = abs( len(seq_lst[-1].strip()) - len("ATC") - 3)              # tail with 3 more nucs
        min_deviation = min( [deviation0, deviation1, deviation2] )
        if deviation1 == min_deviation:
            allele -= 1
        elif deviation2 == min_deviation:
            allele -= 2
        else:
            allele -= 0
    

    if locus in ['D2S441', 'D4S2408', 'D7S820', 'D2S1338', 'D18S51', 'D9S1122', 'vWA']:  # length-based
        seq_lst       = editted_allele.split(' ')[1:-1]
        div = sum([len(i) for i in seq_lst if '[' not in i]) // 4
        mod = sum([len(i) for i in seq_lst if '[' not in i]) % 4
        allele += div + mod / 10
    
    if locus in ['DYS438']:
        seq_lst       = editted_allele.split(' ')[1:-1]
        div = sum([len(i) for i in seq_lst if '[' not in i]) // 5
        mod = sum([len(i) for i in seq_lst if '[' not in i]) % 5
        allele += div + mod / 10

    if locus == 'DYS390':
        allele -= 2

    if locus == 'DYS576':
        allele -= 1   # pat AA[AAAG]1 actually is AAA (in repeat mask) | AAAG, AA might for basecalling-error

    if locus == "D5S818":
        seq_lst    = editted_allele.split(' ')

        if 'CTCT' not in seq_lst[0].strip():
            allele -= 1
        # if abs( len(seq_lst[0].strip()) - len("ATAC") ) < abs( len(seq_lst[0].strip()) - len("ATACCTCT") ):
            # allele -= 1
    
    if locus == 'PentaD':
        if len(unit_cnt_lst) > 2:
            allele = unit_cnt_lst[1] + unit_cnt_lst[2]
        else:
            allele = sum(unit_cnt_lst)

    if  locus.upper() == "DYS612":
        cnts = findall("\d+", editted_allele)
        allele = int(cnts[0]) + int(cnts[-1])

    if locus.upper() == "DYS448":
        cnts = findall("\d+", editted_allele)
        if len(cnts) == 1:
            allele = ""
        else:
            allele = int(cnts[0]) + int(cnts[-1])

    if locus.upper() == "DYS643":
        cnts = findall("\d+", editted_allele)
        allele = int(cnts[0])
    
    if locus.upper() == 'D20S482':
        seq_lst = editted_allele.split(' ')
        if 'AGCT' in seq_lst[-1]:
            allele += 1



    if locus.upper() == 'D19S433':
        seq_lst = editted_allele.split(' ')[1:-1]
        for ind, unit in enumerate(seq_lst):
            if '[' in unit:
                continue
            if 'TT' == unit:
                allele = allele - 1 + 0.2

    if locus.upper() == 'D2S441':
        if editted_allele.split(' ')[0][-1] == 'A':
            allele = allele + 0.1
        else:
            seq_lst = editted_allele.split(' ')[1:-1]
            for ind, unit in enumerate(seq_lst):
                if '[' in unit:
                    continue
                if len(unit) == 4:
                    allele += 1
                else:
                    allele += 0.1 * len(unit)
                    
    
    if locus.upper() == 'D6S1043':
        seq_lst = editted_allele.split(' ')[1:-1]
        for ind, unit in enumerate(seq_lst):
            if '[' in unit:
                continue
            allele = allele + len(unit) * 0.1

    if locus.upper() == 'D1S1656':
        seq_lst = editted_allele.split(' ')[1:-1]
        for ind, unit in enumerate(seq_lst):
            if '[' in unit:
                continue
            if 'TCA' == unit:
                allele = allele + 0.3
    
    if locus.upper() == 'TH01':
        seq_lst = editted_allele.split(' ')[1:-1]
        for ind, unit in enumerate(seq_lst):
            if '[' in unit:
                continue
            if 'ATG' == unit:
                allele = allele + 0.3

    if locus.upper() == 'D21S11':
        seq_lst = editted_allele.split(' ')[1:-1]

        digits = -11
        for ind, unit in enumerate(seq_lst):
            if '[' in unit:
                continue
            digits += len(unit)
        if digits > 0:
            allele += digits/10



    ## digital genotype
    if locus.upper() in ['FGA', 'CSF1PO']:
        pattern = r'[\[ATCG\]]+'
        matches = findall(pattern, editted_allele)

        tot_incomp_len = 0

        for ind, unit in enumerate(matches):
            if (ind == 0 or ind == len(matches) - 1):
                continue
            if '[' in unit:
                continue
            tot_incomp_len += len(unit)
        allele += tot_incomp_len // 4 +  round(tot_incomp_len % 4 * 0.1, 1)
    return allele





def allele_corrector(bracket_allele, ordered_pat_deno_list, counted_pat_dict):
    
    for unit in counted_pat_dict:
        unit_len = len(unit)
        break

    for pat in ordered_pat_deno_list:
        if 'x' in pat:
            return bracket_allele

    bracket_allele_lst = bracket_allele.split()
    bracket_allele_upd = []
    for ind, unit in enumerate(bracket_allele_lst):
        if '[' in unit or ind ==0 or ind == len(bracket_allele_lst) - 1:
            bracket_allele_upd += [unit]
        elif len(unit) == unit_len:
            bracket_allele_upd += [ "[%s]%d" % (unit, 1) ]
        elif len(unit) < unit_len:
            bracket_allele_upd += [ "[%s]0.%d" % (unit, len(unit)) ]
        else:
            bracket_allele_upd += [ unit ]

    return " ".join(bracket_allele_upd)