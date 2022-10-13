#!/usr/bin/env python3

import pandas as pd

def f_shift(row):
    if row['Pos'] <= 248:
        val = row['Pos']+25
    elif row['Pos'] >= 250:
        val = row['Pos']+24
    return val

def main(fname_all_mutations, fname_all_mutations_shifted):

    df_all_muts = pd.read_csv(fname_all_mutations)

    sample = str(fname_all_mutations).split("/")[-1]

    if sample != "parental_stock_ref_EBref":
        df_all_muts = df_all_muts[df_all_muts['Pos']!=249]
        df_all_muts['Pos'] = df_all_muts.apply(f_shift, axis=1)

    df_all_muts.to_csv(fname_all_mutations_shifted)

if __name__ == "__main__":
    main(
        snakemake.input.fname_all_mutations,
        snakemake.output.fname_all_mutations_shifted,
    )
