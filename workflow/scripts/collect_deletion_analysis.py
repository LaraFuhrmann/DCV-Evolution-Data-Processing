#!/usr/bin/env python3
import pandas as pd


def main(in_fnames, out_fname):
    merged_div_csv = pd.concat(
        [pd.read_csv(path_div) for path_div in in_fnames]
    )
    merged_div_csv.to_csv(out_fname)

if __name__ == "__main__":
    main(
        snakemake.input.fnames_rdrp + snakemake.input.fnames_utr,
        snakemake.output.fname,
    )
