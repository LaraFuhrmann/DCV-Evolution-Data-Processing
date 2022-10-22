#!/usr/bin/env python3
"""
Script aggregating co-occurring mutations from all samples.
"""
import pandas as pd
import get_cooccurring_mutations

def main(fname_reference, dnames_shorah, fname_cooccurring_mutations_csv):

    tmp = []

    for shorah_dir in dnames_shorah:

        shorah_dir = str(shorah_dir).split('snvs.vcf')[0]+"REGION_1/"

        sample = shorah_dir.split('/')[-8]
        patient = shorah_dir.split('/')[-6]
        date = shorah_dir.split('/')[-5]

        get_cooccurring_mutations.main(fname_reference, shorah_dir)

        df_tmp = pd.read_csv(shorah_dir+'cooccurring_mutations.csv')
        df_tmp['sample'] = sample
        df_tmp['patient'] = patient
        df_tmp['date'] = date

        tmp.append(df_tmp)

    merged_csv = pd.concat( tmp )
    merged_csv.to_csv(fname_cooccurring_mutations_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fname_reference,
        snakemake.input.dnames_shorah,
        snakemake.output.fname_cooccurring_mutations,
    )
