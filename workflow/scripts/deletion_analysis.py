#!/usr/bin/env python3
"""
Analysis of the occurrence of deletions among the mutation calls.
The purpose of this script is to:
- find locations/regions where a lot of deletions occurr
- find biases of the mutation caller
- understand if Shorah-deletion-problem is appparent in this dataset (meanig very
  low frequent deletions that are long - due to misalignment. This occurs mainly
  with bwa as alginer.)

ShoRAH - Deletion calling:
deletions are reported position-wise.

HIV-1 LTE deletion Analysis
- we expect large deletions in the NEF-region as this part is not used for replication
- Karin: large deletions at the two ends
- we should not have large deletions in other genes (gag, pol) as the virus is still replicating and those are involved in the replication process


Output:
- file with all the deletions listed length and frequences and co.

"""

import pandas as pd
import numpy as np


def ranges(nums):
    """
    auxiliary function for list_frameshift_dels().
    Input is a list of numbers
    return ranges that are covered by those numbers,
    e.g. [1,2,3,10]--> [(1,3),(10,10)]
    """
    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


def len_del(item_range):
    """
    auxiliary function for list_frameshift_dels().
    computing the lenght of item_range,
    """
    return item_range[1] - item_range[0] + 1


def df_chromosome(row):
    return row["Chromosome"][0]


def df_ref(row):
    return "".join(row["Ref"])


def df_var(row):
    return "".join(row["Var"])


def df_length(row):
    return len(row["Var"])


def df_Average_Fvar(row):
    return sum([float(x) for x in row["Fvar"]]) / len(row["Fvar"])


def df_Average_Rvar(row):
    return sum([float(x) for x in row["Rvar"]]) / len(row["Rvar"])


def df_Average_Ftot(row):
    return sum([float(x) for x in row["Ftot"]]) / len(row["Ftot"])


def df_Average_Rtot(row):
    return sum([float(x) for x in row["Rtot"]]) / len(row["Rtot"])


def df_Average_Frequency(row):
    return (row["Average_Fvar"] + row["Average_Rvar"]) / (row["Average_Ftot"] + row["Average_Rtot"])


def parse_deletions_from_csv(fname_snv_csv):
    # Parse snv.csv
    df_snv = pd.read_csv(fname_snv_csv)

    # Filter for deletions
    df_snv = df_snv[df_snv["Var"] == "-"]

    # get length and start positions of deletions
    # and merge the rows.
    deletion_ranges = ranges(df_snv["Pos"])
    deletion_position_length = []
    columns_df_snv = df_snv.columns

    for deletion_range in deletion_ranges:
        deletion_position_length.append(
            [deletion_range[0], deletion_range[1], len_del(deletion_range)]
        )

    df_snv["Start_Position"] = -1

    for row_iter, row in df_snv.iterrows():
        for deletion_range in deletion_ranges:
            if row["Pos"] in range(deletion_range[0], deletion_range[1] + 1):
                df_snv.at[row_iter, "start_position"] = deletion_range[0]

    df_snv = (
        df_snv.groupby("start_position")[columns_df_snv]
        .agg(lambda x: list(x))
        .reset_index()
    )

    # formatting and post-Processing
    df_snv["Chromosome"] = df_snv.apply(df_chromosome, axis=1)
    df_snv["Ref"] = df_snv.apply(df_ref, axis=1)
    df_snv["Var"] = df_snv.apply(df_var, axis=1)
    df_snv["Length"] = df_snv.apply(df_length, axis=1)
    df_snv["Average_Fvar"] = df_snv.apply(df_Average_Fvar, axis=1)
    df_snv["Average_Rvar"] = df_snv.apply(df_Average_Rvar, axis=1)
    df_snv["Average_Ftot"] = df_snv.apply(df_Average_Ftot, axis=1)
    df_snv["Average_Rtot"] = df_snv.apply(df_Average_Rtot, axis=1)

    df_snv["Average_Frequency"] = df_snv.apply(df_Average_Frequency, axis=1)

    df_snv = df_snv.drop("Pos", axis=1)

    return df_snv


def main(fname_snv_csv, fname_output):

    cell_line = fname_snv_csv.split("/variants")[0].split("/")[-2]
    passage = fname_snv_csv.split("/variants")[0].split("/")[-1]

    # Parse snv.csv
    df_deletions = parse_deletions_from_csv(fname_snv_csv)
    df_deletions["cell_line"] = cell_line
    df_deletions["passage"] = passage
    df_deletions.to_csv(fname_output)


    # check for long, very low frequent deletions

    # plot with deletions marked along the genome

if __name__ == "__main__":
    main(
        snakemake.input.fname_snv_csv,
        snakemake.output.fname_out,
    )
