#!/usr/bin/env python3
import subprocess
from fuc import pyvcf


def transform_to_EB_space(in_vcf, out_vcf, sample):
    """
    adapt positions of mutation calls such that all is wrt the EB-reference space
    parental_stock.conensus.fasta.
    We are adjusting the parental_stock.consensus.bcftools.fasta sequence such
    that it excludes the UTRs.
    This means we exclude positions [0,24] and [9102, 9263] --
    those positions all have less than 9999 coverage.
    --> to map back we need to add 25 to the POS column.

    ---
    There is an insertion in the parental_stock.consensus.bcftools.fasta at
    position 273, we have T --> TG.
    To move back to the EB ref space, we do the following:

    248 --> +25 = 273 [corresponds to T in EB ref]
    249 --> delete [this is where the insertion "G" is]
    250 --> +24 = 274 [corresponds to A in EB ref]
    """
    vf = pyvcf.VcfFrame.from_file(in_vcf)
    vf.df["CHROM"] = "AF014388"
    # exclude deletions since vcf-annotator does not know how to treat them in this form.
    vf.df = vf.df[vf.df["ALT"] != "-"]

    if sample != "parental_stock_ref_EBref":
        # exclude mutations at position 249
        vf.df = vf.df[vf.df["POS"] != 249]
        vf.df["POS"] = vf.df.apply(f_shift, axis=1)

    vf.to_file(out_vcf)


def f_shift(row):
    if row["POS"] <= 248:
        val = row["POS"] + 25
    elif row["POS"] >= 250:
        val = row["POS"] + 24
    return val


def add_vcf_DP_and_AF(vcf_in, vcf_out):
    """
    # alter vcf file such that it matches the necessary input of SNPGenie
    # https://github.com/chasewnelson/SNPGenie#snpgenie-input
    """
    vf = pyvcf.VcfFrame.from_file(vcf_in)

    Rtot = vf.df["INFO"].str.split("Rtot=").str[1].str.split(";").str[0]
    Ftot = vf.df["INFO"].str.split("Ftot=").str[1].str.split(";").str[0]
    DP = Rtot.astype("int") + Ftot.astype("int")

    Rvar = vf.df["INFO"].str.split("Rvar=").str[1].str.split(";").str[0]
    Fvar = vf.df["INFO"].str.split("Fvar=").str[1].str.split(";").str[0]
    AF = (Rvar.astype("int") + Fvar.astype("int")) / DP

    vf.df["INFO"] = vf.df["INFO"] + ";DP=" + DP.astype("str")
    vf.df["INFO"] = vf.df["INFO"] + ";AF=" + AF.astype("str")

    vf.to_file(vcf_out)


def run_snpgenie(fname_reference, in_vcf, gtffile_CDS_annotations, dname_work):
    """
    WITHIN-POOL ANALYSIS. Use snpgenie.pl, the original SNPGenie.
    Analyzes within-sample πN/πS from pooled NGS SNP data.

    fname_reference: fasta of reference sequence wrt which the SNVs where called
    """
    subprocess.run(
        [
            "snpgenie.pl",
            "--vcfformat=2",
            "--snpreport=" + str(in_vcf),
            "--fastafile=" + str(fname_reference),
            "--gtffile=" + str(gtffile_CDS_annotations),
            #"--minfreq=0.001",
            "--outdir="+str(dname_work)
        ],
        check=True,
        #cwd=dname_work,
    )


def main(fname_reference, fname_snv_in, gtffile_CDS_annotations, dname_work):

    # map vcf file to EB ref
    sample = str(fname_snv_in).split("/variants")[0].split("/")[-4]
    fname_snv_temp = str(fname_snv_in).split(".vcf")[0] + ".temp.snpgenie.vcf"
    transform_to_EB_space(fname_snv_in, fname_snv_temp, sample)

    # add DP4 to vcf file
    add_vcf_DP_and_AF(fname_snv_temp, fname_snv_temp)

    from pathlib import Path
    import shutil

    test = Path(dname_work)
    if test.exists() and test.is_dir():
        shutil.rmtree(test)

    # run snpgenie
    run_snpgenie(fname_reference, fname_snv_temp, gtffile_CDS_annotations, dname_work)


if __name__ == "__main__":
    main(
        snakemake.input.fname_reference,
        snakemake.input.fname_snvs_vcf,
        snakemake.input.fname_gtffile,
        snakemake.output.dname_work
    )
