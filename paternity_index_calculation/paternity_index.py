"""
This is the code for Paternity Index calculation

"""
import pandas as pd 
from collections import defaultdict

class Paternity_Index_Cal:
    def __init__(self, freq_csv_path):
        self.freq_csv_path = freq_csv_path
        self.allele_freq_dict = self.get_allele_freq_dict()

    def single_paternity_cal(self, child_alleles, pater_alleles, locus):
        child_allele_1, child_allele_2 = child_alleles
        pater_allele_1, pater_allele_2 = pater_alleles
        locus_allele_freq_dict         = self.allele_freq_dict[locus]

        if child_allele_1 == child_allele_2 and pater_allele_1 == pater_allele_2:

            if child_allele_1 == pater_allele_1:
                return 1/locus_allele_freq_dict[child_allele_1] # q, q, 1/q
            else:
                return 0  # q,q

        if child_allele_1 == child_allele_2 and pater_allele_1 != pater_allele_2:
            if child_allele_1 == pater_allele_1 or child_allele_1 == pater_allele_2:
                return 1/(2*locus_allele_freq_dict[child_allele_1]) # q,qr
            else:
                return 0 # q, pr
        
        if child_allele_1 != child_allele_2 and pater_allele_1 == pater_allele_2:
            if child_allele_1 == pater_allele_1:
                return 1/(2*locus_allele_freq_dict[child_allele_1]) # q,qr
            elif child_allele_2 == pater_allele_1:
                return 1/(2*locus_allele_freq_dict[child_allele_2])  # pq, q
            else:
                return 0 # pq, r

        if child_allele_1 != child_allele_2 and pater_allele_1 != pater_allele_2:
            if (child_allele_1 == pater_allele_1 and child_allele_2 == pater_allele_2) or (child_allele_1 == pater_allele_2 and child_allele_2 == pater_allele_1):
                return (locus_allele_freq_dict[child_allele_1] + locus_allele_freq_dict[child_allele_2]) / (4 * locus_allele_freq_dict[child_allele_1] * locus_allele_freq_dict[child_allele_2])  # pq, pq, (p+q) / 4pq	
            
            elif  child_allele_1 == pater_allele_1 or child_allele_1 == pater_allele_2:
                return 1/(4*locus_allele_freq_dict[child_allele_1]) # pq, qr, 1/4q
            
            elif child_allele_2 == pater_allele_1 or child_allele_2 == pater_allele_2 :
                return 1/(4*locus_allele_freq_dict[child_allele_2]) # pq, qr
            else:
                return 0 # pq, rs
    
    def get_allele_freq_dict(self):
        
        freq_dict = defaultdict(lambda: defaultdict(int))

        with open(self.freq_csv_path, 'r') as f:
            for line_cnt, line in enumerate(f):
                if line_cnt == 0:
                    continue
                locus, allele, freq = line.strip().split(',')
                freq_dict[locus][float(allele)] = float(freq)
        return freq_dict



class LR_Index_Cal:
    def __init__(self, freq_csv_path):
        self.freq_csv_path = freq_csv_path
        self.allele_freq_dict = self.get_allele_freq_dict()

    def individual_identification_cal(self, alleles, locus):
        if pd.isna(alleles):
            return pd.NA

        allele_1, allele_2 = float(alleles.split(',')[0]), float(alleles.split(',')[1])
        locus_allele_freq_dict         = self.allele_freq_dict[locus]

        if allele_1 == allele_2:
            return locus_allele_freq_dict[allele_1] ** 2

        else:
            return 2 * locus_allele_freq_dict[allele_1] * locus_allele_freq_dict[allele_2]
        
    
    def get_allele_freq_dict(self):
        
        freq_dict = defaultdict(lambda: defaultdict(int))

        with open(self.freq_csv_path, 'r') as f:
            for line_cnt, line in enumerate(f):
                if line_cnt == 0:
                    continue
                locus, allele, freq = line.strip().split(',')
                freq_dict[locus][float(allele)] = float(freq)
        return freq_dict

if __name__ == '__main__':
    PIC = Paternity_Index_Cal('CHN_STR_allele_freqs.csv')
    print( PIC.single_paternity_cal([13.0, 10.0], [12.0, 10.0], 'CSF1PO') )

    PIC = LR_Index_Cal('CHN_STR_allele_freqs.csv')
    print( PIC.individual_identification_cal([13.0, 10.0], 'CSF1PO') )
    
