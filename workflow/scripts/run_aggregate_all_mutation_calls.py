import os
import pandas as pd


base_path='/cluster/work/bewi/members/lfuhrmann/HIV-LTE/HIV-LTE-NGS-data-Experiment3/samples/'

cell_lines = ['MT2_1','MT2_2','MT4_1','MT4_2']
viral_passage = ['VP'+str(i) for i in range(10,340,10)]

file_path='/variants/SNVs/'
in_fnames_snv_vcf =[]

for MT in cell_lines:
    for passage in viral_passage:
        fname_snv_in= base_path+MT+'/'+passage+file_path+'snvs.csv'
        if os.path.isfile(fname_snv_in)==True:
            in_fnames_snv_vcf.append(fname_snv_in)
        else:
            print('Not found: ',MT,'/',passage)
            print('path: ',fname_snv_in)

import aggregate_mutation_calls

aggregate_mutation_calls.main(in_fnames_snv_vcf, '/cluster/work/bewi/members/lfuhrmann/HIV-LTE/HIV-LTE-NGS-data-Experiment3/all_mutations.csv' )
