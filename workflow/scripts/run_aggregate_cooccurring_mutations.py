import os
import pandas as pd

fname_reference='/cluster/work/bewi/members/lfuhrmann/HIV-LTE/HIV-LTE-NGS-data-Experiment3/references/ancestor_consensus.fasta'
base_path='/cluster/work/bewi/members/lfuhrmann/HIV-LTE/HIV-LTE-NGS-data-Experiment3/samples/'

cell_lines = ['MT2_1','MT2_2','MT4_1','MT4_2']
viral_passage = ['VP'+str(i) for i in range(10,340,10)]

dnames_shorah =[]

for MT in cell_lines:
    for passage in viral_passage:
        fname_snv_in= base_path+MT+'/'+passage+'/variants/SNVs/REGION_1/'
        if os.path.isdir(fname_snv_in)==True:
            dnames_shorah.append(fname_snv_in)
        else:
            print('Not found: ',MT,'/',passage)
            print('path: ',fname_snv_in)

import aggregate_cooccurring_mutations

fname_cooccurring_mutations_csv = '/cluster/work/bewi/members/lfuhrmann/HIV-LTE/HIV-LTE-NGS-data-Experiment3/all_cooccurring_mutations.csv'
aggregate_cooccurring_mutations.main(fname_reference, dnames_shorah, fname_cooccurring_mutations_csv)
