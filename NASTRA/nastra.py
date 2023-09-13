from pathlib            import Path
from collections        import Counter, defaultdict
from libs.config_loader import locus_info_loader, pattern_loader, bam_loader
from libs.pairwise_alignment import bioalign, read_trim, read_group, search_tree
from libs.reads_genotyper import reads_genotyper
from libs.genotype_merge import merged_func

import pandas as pd
from tqdm import tqdm
import os, sys, argparse


def main():
    args = get_args()
    args.func(args)


def get_args():
    script_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(description='NanoSTR: a STR genotyping tool for forensic STR')
    subparsers = parser.add_subparsers(help='command: {call, }')
    subparsers.required = True

    ## call
    main_parser = subparsers.add_parser('call',   help='output calling result')
    main_parser.add_argument('-b', '--bam',       help='the path of bam file', required = True, type = str)
    main_parser.add_argument('-f', '--factsheet', help='the path of STR factsheet', default = "%s/cfgs/repeat_structure.pat" % script_dir, required = False, type = str)
    main_parser.add_argument('-p', '--panel',     help='the path of the loci panel file', default = "%s/cfgs/panel_forenseq.csv" % script_dir, required=False, type = str)
    main_parser.add_argument('-c', '--config',    help='the config file of SN cutoff', default = "%s/cfgs/threshold.cfg" % script_dir, required=False, type=str)
    main_parser.add_argument('-o', '--output',    help='the path of output', required = True, type = str)
    main_parser.add_argument('--samtools',        help='the PATH to samtools', default = "samtools", required = False, type = str)
    main_parser.add_argument('--sncutoff',        help='the value of sn cutoff', default = 25, required = False, type = int)
    main_parser.set_defaults(func=calling_func)

    ## run
    args = parser.parse_args()
    return args


def genotype_calling(nanostr_lst, sn_cutoff, snr_cutoff):
    if len(nanostr_lst) == 0:
        return 'Fail: NoCalling'
    elif len(nanostr_lst) == 4:
        nst1, nst_sn1, nst_snr1, nst_seq1 = nanostr_lst
        if nst_sn1 <= sn_cutoff:
            return 'Fail: Interpretation'
        nst2 = nst1
        nst_sn2  = 0
        nst_seq2 = nst_seq1
    else:
        nst1, nst_sn1, nst_snr1, nst_seq1, nst2, nst_sn2, nst_snr2, nst_seq2 = nanostr_lst
        if nst_sn1 < sn_cutoff and nst_sn2 < sn_cutoff:
            return 'Fail: Interpretation'
        if nst_sn1 >= sn_cutoff and nst_sn2 < sn_cutoff:
            return 'Fail: Imbalance'
        if nst_snr2 < snr_cutoff:
            nst2     = nst1
            nst_sn2  = 0
            nst_seq2 = nst_seq1
    return nst1, nst2, nst_sn1, nst_sn2, round(nst_sn2/nst_sn1, 3), nst_seq1, nst_seq2, 


def calling_func(args):
    bamfile         = Path(args.bam)
    cfgpath         = Path(args.config)
    fact_sheet_path = Path(args.factsheet)
    locus_info_file = Path(args.panel)
    samtools        = args.samtools
    sn_cutoff       = args.sncutoff
    outpath         = Path(args.output)
    logout          = outpath.parent / (outpath.stem + '.log')

    barcode         = bamfile.stem
    locus_info_dict = locus_info_loader(locus_info_file)
    pairwise        = bioalign(75, -90, 75, 10)
    cluster_func    = read_group(aln_func = pairwise, num_candidates = 6, max_num_groups = 200, cluster_score_thres = 1)
    
    results = []
    for ind, locus in tqdm(enumerate(locus_info_dict), desc = '[Processing]'):
    # for ind, locus in enumerate(locus_info_dict):
        counted_pat_dict, ordered_pat_deno_list, ordered_pat_list = pattern_loader(fact_sheet_path, locus)
        region, flank_len, prefix, suffix = locus_info_dict[locus]
        
        search_func  = search_tree(counted_pat_dict)
        trim_func    = read_trim(pairwise, prefix, suffix, reserved_flanking_len=3)

        raw_read_ids, raw_reads = bam_loader(samtools, bamfile, region, depth=-1)

        if len(raw_reads) < 10:
            result = [[barcode, locus, '<None>', '<None>', len(raw_reads), '<LowCov>']]
            clipped_dat = pd.DataFrame(result, columns=['barcode', 'locus', 'seq', 'genotype', 'sn', 'sn_ratio'])
            results.append(clipped_dat)
            continue

        trimmed_reads           = [trim_func.trim(read) for read in raw_reads]
        counter_dct             = Counter(trimmed_reads)
        
        cluster_alleles         = cluster_func.cluster(counter_dct)
        brac_alleles            = [(search_func.search(seq), cnt) for seq, cnt in cluster_alleles]
        result = []
        for brac_allele, brac_num in brac_alleles:
            if brac_allele is None or brac_allele == "" or "[" not in brac_allele: # brac_num -> sn <=10 is invalid
                result.append( [barcode, locus, '<None>', '<None>', brac_num] )
            else:
                editted_brac_allele, repeat_count = reads_genotyper(locus, brac_allele)
                result.append( [barcode, locus, brac_allele, repeat_count, brac_num] )
        
        result = merged_func(result, pairwise)
        
        clipped_dat = pd.DataFrame(result, columns=['barcode', 'locus', 'seq', 'genotype', 'sn']).sort_values('sn', ascending=False)
        
        if len(clipped_dat) == 0:
            continue
        else:
            clipped_dat['sn_ratio'] = round( clipped_dat['sn'] / max( clipped_dat['sn'] ), 5 )
            results.append(clipped_dat)

    merged_dat = pd.concat(results, axis=0)
    logout.parent.mkdir(exist_ok=True, parents=True)
    merged_dat.to_csv(logout, index=None)

    ## calling
    best_threshold_df = pd.read_csv(cfgpath)
    genotype_dict     = defaultdict(list)

    for row in merged_dat.values:
        barcode, locus, seq, genotype, sn,sn_ratio = row
        if len( genotype_dict['%s_%s' % (barcode, locus)] ) > 6:
            continue
        if genotype == '<None>' or genotype == '':
            continue        
        genotype_dict['%s_%s' % (barcode, locus)] += [float(genotype), int(sn), float(sn_ratio), seq]

    outfile = open(outpath, 'w')
    outfile.write('barcode,locus,qc_info,allele1,allele2,sn1,sn2,snr,seq1,seq2\n')
    for key in genotype_dict:
        sample_name, locus = key.split('_')
        snr_cutoff     = best_threshold_df.loc[best_threshold_df.locus == locus, 'cov_%d' % sn_cutoff].values[0]
        nanostr_out    = genotype_calling(genotype_dict[key], sn_cutoff, snr_cutoff)

        if 'Fail' in nanostr_out:
            outfile.write( f"{sample_name},{locus},{nanostr_out},,,,,,,\n" )
        else:
            outfile.write( f"{sample_name},{locus},pass,{','.join([str(i) for i in nanostr_out])}\n" )
    outfile.close()




if __name__ == '__main__':
    # print(os.path.split(os.path.realpath(__file__))[0])
    main()
