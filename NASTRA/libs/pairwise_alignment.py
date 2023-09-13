import parasail
from collections import defaultdict
import re
"""
1. trim with 3bp flanking reserved
2. cluster

"""

class bioalign:
    def __init__(self, match, mismatch, gap_open, gap_ext):
        self.user_matrix = parasail.matrix_create("ACGT", match, mismatch)
        self.gap_open    = gap_open
        self.gap_ext     = gap_ext

    def align(self, refer, query):
        # aln = parasail.sg_trace(query, refer, self.gap_open, self.gap_ext, self.user_matrix)

        ## Smith-Waterman local alignment
        aln = parasail.sw_trace(query, refer, self.gap_open, self.gap_ext, self.user_matrix)
        return aln

    def align_nw(self, refer, query):
        ## Needleman-Wunsch global alignment
        aln = parasail.nw_trace(query, refer, self.gap_open, self.gap_ext, self.user_matrix)
        return aln
        

class read_trim:
    def __init__(self, aln_func, pre_ref, suf_ref, reserved_flanking_len = 3):
        self.len = reserved_flanking_len
        self.pairwise = aln_func
        self.pre_ref  = pre_ref
        self.suf_ref  = suf_ref

    def trim(self, query):
        trimmed_1 = self.trim_head(query)
        trimmed_2 = self.trim_tail(trimmed_1)
        return trimmed_2

    def trim_head(self, query):
        aln = self.pairwise.align(self.pre_ref, query)
        start_point = max(aln.end_query-self.len+1, 0)
        trimmed_read = query[start_point:]
        return trimmed_read

    def trim_tail(self, query):
        aln = self.pairwise.align(self.suf_ref[::-1], query[::-1])
        start_point = max(aln.end_query-self.len+1, 0)
        trimmed_read  = query[::-1][(start_point):]
        return trimmed_read[::-1]




class read_group:
    def __init__(self, aln_func, num_candidates = 3, max_num_groups = 200, cluster_score_thres = 2):
        self.pairwise = aln_func
        self.max_num_groups  = max_num_groups
        self.num_candidates  = num_candidates
        self.cluster_score_thres = cluster_score_thres

    def cluster(self, counter_dct):
        part_group = counter_dct.most_common(self.max_num_groups)
        allele_dct = self.allele_init(part_group)
        for allele, supnum in part_group[1:]:
            group_name, group_count = self.is_new_allele(allele_dct, allele, supnum)
            if group_name is None:
                continue
            allele_dct[group_name] += group_count
        
        return [(k, v) for k, v in allele_dct.items()]

    def allele_init(self, part_group):
        allele_dct         = defaultdict(int) 
        allele, supnum     = part_group[0]
        allele_dct[allele] = supnum
        return allele_dct

    def is_new_allele(self, allele_dct, allele, supnum):

        if len(allele_dct) < self.num_candidates:
            for candidate in allele_dct:
                aln       = self.pairwise.align_nw(candidate, allele)
                cigar_str = aln.cigar.decode.decode('UTF-8')
                sym_str   = ''.join((re.findall('[DIX=]+', cigar_str)))
                ## same allele
                # print(candidate)
                # print(allele)
                # print(cigar_str, self.get_cigar_score(cigar_str))
                # print()
                if self.get_cigar_score(cigar_str) <= 1: # one gap in core region
                    if sym_str.startswith('I') or sym_str.endswith('I') or sym_str.startswith('D') or sym_str.endswith('D'): # require the I/D must in head or tail
                        return candidate, supnum
                    
                    elif 'X' in sym_str:
                        return allele, supnum

                    else:  # ignore if
                        return candidate, supnum ## internal indel also assigned as same allele
            ## new allele
            return allele, supnum
        
        else:
            best_class = None
            best_score = 1e6
        
            for candidate in allele_dct:
                aln         = self.pairwise.align_nw(candidate, allele)
                cigar_str   = aln.cigar.decode.decode('UTF-8')
                cigar_score = self.get_cigar_score(cigar_str)

                if cigar_score < best_score:
                    best_score = cigar_score
                    best_class = candidate

            if best_score > self.cluster_score_thres:
                return None, supnum
            else:
                return best_class, supnum
    
    def get_cigar_score(self, cigar_str):
        numIs = re.findall(r'(\d+)I', cigar_str)
        numDs = re.findall(r'(\d+)D', cigar_str)
        numXs = re.findall(r'(\d+)X', cigar_str)

        score = 0
        for i in numIs + numDs + numXs:
            score += 10**( int(i) - 1 )
        return score


class search_tree:
    def __init__(self, pattern_dct):
        self.pattern_dct = pattern_dct

    def search(self, query):
        full_seq_set        = {query}    
        stillContinue       = 1

        while stillContinue:
            full_seq_set_update = self.searchTree(full_seq_set, self.pattern_dct)
            stillContinue       = self.check_stop(full_seq_set, full_seq_set_update)
            full_seq_set        = full_seq_set_update
            # print(len(full_seq_set))
            
        if len(full_seq_set) > 1e5:
            return None 
        best_brac_seq = self.find_best_bracket(full_seq_set)
        for key in self.pattern_dct:
            symbol = self.pattern_dct[key]
            best_brac_seq = best_brac_seq.replace(symbol, key)
        return best_brac_seq.strip()

    def searchTree(self, seq_set, pattern_dict):
        new_seq_lst = []
        
        for seq in seq_set:
            no_update = 1
            for pattern in pattern_dict:
                index_lst = [list(i.span()) for i in re.finditer(pattern, seq)]

                if index_lst == []:
                    continue
                else:
                    pattern_span_lst, unit_count_lst = self.combine_interval(index_lst)

                    new_seqs = [ seq[:pattern_span[0]] \
                                + ' [%s]%d ' % ( pattern_dict[pattern], unit_count_lst[ind] ) \
                                + seq[pattern_span[1]:] for ind, pattern_span in enumerate(pattern_span_lst) ] 

                    new_seq_lst += new_seqs
                    no_update = 0
                    
            if no_update:
                new_seq_lst.append(seq)
        return set(new_seq_lst)

    def combine_interval(self, index_lst):
        unit_count_list  = [1]
        pattern_span_lst = [index_lst[0]]
        
        
        for span_index in index_lst[1:]:
            
            if pattern_span_lst[-1][-1] == span_index[0]:
                pattern_span_lst[-1][-1] = span_index[1]
                unit_count_list[-1]     += 1
            else:
                pattern_span_lst.append(span_index)
                unit_count_list.append(1)
            
        return pattern_span_lst, unit_count_list

    def find_best_bracket(self, full_seq_set):
        
        best_score = 1e6
        best_brac  = None
        for seq in full_seq_set:
            pattern = r'[\[ATCG\]]+'
            matches = re.findall(pattern, seq)
            
            extra_bases = 0
            extra_parts = 0

            for ind, unit in enumerate(matches):
                if '[' in unit:
                    continue
                if (ind == 0 or ind == len(matches) - 1) and len(unit) <= 3:
                    continue
                extra_bases += len(unit)
                extra_parts += 1

            score = extra_parts + 10 * extra_bases
            if score <= best_score:
                best_brac  = seq
                best_score = score
                
        return best_brac

    def check_stop(self, full_seq_set, full_seq_set_update):

        if full_seq_set == full_seq_set_update:
            return 0
        elif len(full_seq_set_update) > 1e5:
            return 0
        else:
            return 1

## trim

# from pathlib import Path
# from config_loader import locus_info_loader, pattern_loader, bam_loader

# bamfile = Path('/data0/Projects/Project_nanoSTR/dataset/alignment/Forenseq_fast/pool1/barcode01.bam')
# fact_sheet_path = '/data0/Projects/Project_nanoSTR/nanoSTR/cfgs/vaTable_revised.txt'
# locus_info_file = '/data0/Projects/Project_nanoSTR/nanoSTR/cfgs/panel_forenseq.csv'
# locus   = 'D22S1045'
# samtools = 'samtools'

# barcode = bamfile.stem
# counted_pat_dict, ordered_pat_deno_list, ordered_pat_list = pattern_loader(fact_sheet_path, locus)
# locus_info_dict = locus_info_loader(locus_info_file)
# region, flank_len, prefix, suffix = locus_info_dict[locus]
# raw_read_ids, raw_reads = bam_loader(samtools, bamfile, region, depth=-1)


# print(len(raw_reads))
# print(prefix, suffix)


# pairwise     = bioalign(75, -90, 75, 10)
# trim_func    = read_trim(pairwise, prefix, suffix)
# cluster_func = read_group(pairwise)
# search_func  = search_tree(counted_pat_dict)


# from time import time
# start = time()
# trimmed_reads = [trim_func.trim(read) for read in raw_reads]
# end   = time()
# print('%% <trimmed_reads> %.5f' % (end - start))

# start = time()
# counter_dct   = Counter(trimmed_reads)
# end   = time()
# print('%% <Counter> %.5f' % (end - start))


# start = time()
# cluster_alleles  = cluster_func.cluster(counter_dct)
# end   = time()
# print('%% <Cluster> %.5f' % (end - start))


# start = time()
# brac_alleles = [(search_func.search(seq), cnt) for seq, cnt in cluster_alleles]
# end   = time()
# print('%% <Search> %.5f' % (end - start))


# print(counter_dct.most_common(5))
# print(cluster_alleles)
# print(brac_alleles)








# for ind, (trimed, cnt) in enumerate(dct.most_common(100)):
#     if ind == 0:
#         alleles[trimed] = cnt
#         continue
#     if ind > 3:
#         continue
    
#     for allele in alleles:
#         aln = pairwise.align(allele, trimed)
#         print(ind, "\n")
#         print(aln.cigar.decode.decode('UTF-8'))
#         print(get_cigar_score(aln.cigar.decode.decode('UTF-8')))
#         # print(get_cigar_score(aln.cigar.decode))
#         print(aln.traceback.ref)
#         print(aln.traceback.query)
#         print('\n')
        
















# pre_ref = 'TTAAGGAGAGTGTCACTATA'
# suf_ref = 'AAACACTATATATATATAAC'

# query   = 'TGGGCCAGGCTATGTAGTGAGGACAAGGAGTCCATCTGGGTTAAGTGATAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTATCTATCTATCTAAAACACTATATATATATAACACTATATATATAATACTATATATATATTAAAAAACACCATAACAGAAACTCAGTAGTCATAGTGAAATCAAATGGAATTCTCGGGTGCCAAGGAACTCCA'

# print(trim_func.trim(query, pre_ref, suf_ref))


# print(user_matrix.matrix)

# refer = 'CTGTATTTACAAATACAT'
# query = 'ATCTGTATACATACATTATCTATC'

# refer = refer[::-1]
# query = query[::-1]


# start = time()
# pairwise = bioalign(75, -90, 75, 10)
# aln = pairwise.align(query, suf_refer)
# end   = time()
# print('\n%.5f\n' % (end - start))


# cigar = aln.cigar
# # cigars have seq, len, beg_query, and beg_ref properties
# # the seq property is encoded
# # use decode attribute to return a decoded cigar string

# print(aln.traceback.ref)
# print(aln.traceback.query)
# # print('\n')
# print('query end = {}'.format(aln.end_query))
# print('ref end = {}'.format(aln.end_ref))
# print('query start = {}'.format(cigar.beg_query))
# print('ref start = {}'.format(cigar.beg_ref))

# # print('trimed', query[(aln.end_query+1):])



# from Bio.pairwise2 import align as bioalign
# MATCH_SCORE         = .75
# MISMATCH_SCORE      = -.9
# GAP_OPEN_SCORE      = -.75
# GAP_EXTEND_SCORE    = -.1
# GAP_CHAR            = "-"
# ALIGN_RES_IND       = 0

# start = time()
# [aln_ref, aln_seq, aln_score, aln_start, aln_end] = bioalign.localms(refer.upper(), query.upper(), 
#                                                                         MATCH_SCORE, MISMATCH_SCORE, 
#                                                                         GAP_OPEN_SCORE, GAP_EXTEND_SCORE, 
#                                                                         gap_char = GAP_CHAR)[ALIGN_RES_IND]
# end   = time()


# print('\n%.5f\n' % (end - start))

# print(aln_ref)
# print(aln_seq)
# print(aln_start)
# print(aln_end)
