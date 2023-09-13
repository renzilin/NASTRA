import re

def get_cigar_score(cigar_str):
    numIs = re.findall(r'(\d+)I', cigar_str)
    numDs = re.findall(r'(\d+)D', cigar_str)
    numXs = re.findall(r'(\d+)X', cigar_str)

    score = 0
    for i in numIs + numDs + numXs:
        score += 10**( int(i) - 1 )
    return score

def allele_coverter(brac_allele, trim_flanking = True):
    
    brac_allele_lst = brac_allele.split(' ')

    if trim_flanking:
        if '[' not in brac_allele_lst[0]:
            brac_allele_lst = brac_allele_lst[1:]
        if '[' not in brac_allele_lst[-1]:
            brac_allele_lst = brac_allele_lst[:-1]

    seq = ''
    for brac in brac_allele_lst:
        if '[' not in brac:
            seq += brac
        else:
            unit = re.findall(r'[ATCG]+', brac)[0]
            cnt  = re.findall(r'\d+', brac)[0]
            seq += unit * int(cnt)
    return seq

def merged_func(result, aln_func):
    genotype_dct = {}
    allele_dct   = {}

    result = [each for each in result if '<None>' not in each] 

    for barcode, locus, brac_allele, repeat_count, brac_num in result:
        if repeat_count not in genotype_dct:
            genotype_dct[repeat_count] = brac_allele
            allele_dct[brac_allele]    = [barcode, locus, brac_allele, repeat_count, brac_num]
        
        else:
            aln         = aln_func.align( allele_coverter(genotype_dct[repeat_count]), allele_coverter(brac_allele) )
            # print(locus, genotype_dct[repeat_count], brac_allele)

            cigar_str   = aln.cigar.decode.decode('UTF-8')
            cigar_score = get_cigar_score(cigar_str)
            # print(genotype_dct[repeat_count], brac_allele, cigar_score)
            if cigar_score <= 1:
                allele_dct[ genotype_dct[repeat_count] ][-1] += brac_num
            else:
                allele_dct[brac_allele] = [barcode, locus, brac_allele, repeat_count, brac_num]
    
    output = [ allele_dct[key] for key in allele_dct ] 
    # print(output)
    return output