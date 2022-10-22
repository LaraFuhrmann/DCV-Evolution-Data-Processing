#!/usr/bin/env python3

import pandas as pd

def f_shift(row):
    if row['position'] <= 248:
        val = row['position']+25
    elif row['position'] >= 250:
        val = row['position']+24
    else:
        val = 0
    return val

def main(fname_all_mutations, fname_all_mutations_shifted):

    df_all_muts = pd.read_csv(fname_all_mutations)

    sample = str(df_all_muts['sample'][0])

    if sample != "parental_stock_ref_EBref":
        df_all_muts['position_EB_space'] = df_all_muts.apply(f_shift, axis=1)

    df_all_muts.to_csv(fname_all_mutations_shifted)

if __name__ == "__main__":
    main(
        snakemake.input.fname_all_mutations,
        snakemake.output.fname_all_mutations_shifted,
    )
