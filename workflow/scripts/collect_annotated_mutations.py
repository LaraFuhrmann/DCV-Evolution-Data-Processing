#!/usr/bin/env python3
"""
Collect annotated mutations from all samples.
"""
import pandas as pd
from fuc import pyvcf

def main(fnames_snv_csv, fout_all_mutations_csv):

    tmp = []

    for f_snv_vcf in fnames_snv_csv:

        f_snv_vcf = str(f_snv_vcf)
        df_vcf = pyvcf.VcfFrame.from_file(f_snv_vcf).df

        df_vcf["genotype"] = f_snv_vcf.split("/variants")[0].split("/")[-4]
        df_vcf["replicate"] = f_snv_vcf.split("/variants")[0].split("/")[-2]
        df_vcf["passage"] = f_snv_vcf.split("/variants")[0].split("/")[-1]

        # INFO field to dataframe
        info_strings = '{"' + df_vcf.INFO.str.split(';').str.join('","').str.replace('=','":"').str.replace("\"\",", "") + '"}'
        info_df = pd.json_normalize(info_strings.apply(eval))

        df_tmp = pd.concat([df_vcf, info_df], axis=1)
        tmp.append(df_tmp)

    merged_div_csv = pd.concat(
        tmp
    )
    merged_div_csv.to_csv(fout_all_mutations_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fnames_snv_csv,
        snakemake.output.fname_all_mutations,
    )
