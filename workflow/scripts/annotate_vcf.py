#!/usr/bin/env python3
"""
Annotate mutation calls using https://github.com/rpetit3/vcf-annotator.
Position of mutation calls are with respect to the modified parental_stock_consensus
and therefore must be mapped back to EB reference space.

"""
from fuc import pyvcf
import subprocess

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
    vf.df['CHROM'] = 'AF014388'
    # exclude deletions since vcf-annotator does not know how to treat them in this form.
    vf.df = vf.df[vf.df['ALT']!='-']

    if sample != "parental_stock_ref_EBref":
        # exclude mutations at position 249
        vf.df = vf.df[vf.df['POS']!=249]
        vf.df['POS'] = vf.df.apply(f_shift, axis=1)

    vf.to_file(out_vcf)

def f_shift(row):
    if row['POS'] <= 248:
        val = row['POS']+25
    elif row['POS'] >= 250:
        val = row['POS']+24
    return val

def run_vcf_annotator(in_vcf, path_vcf_annotator, genbank_file, out_vcf):

    subprocess.run(
        [
            "python3",
            str(path_vcf_annotator),
            str(in_vcf),
            str(genbank_file),
            "--output",
            out_vcf,
        ],
        check=True,
    )

def main(fname_snv_in, path_vcf_annotator, fname_genbank_file, fname_snv_out):

    sample = str(fname_snv_out).split("/variants")[0].split("/")[-4]
    fname_snv_temp = str(fname_snv_in).split('.vcf')[0]+'.temp.vcf'

    transform_to_EB_space(fname_snv_in, fname_snv_temp, sample)


    run_vcf_annotator(fname_snv_temp, path_vcf_annotator, fname_genbank_file, fname_snv_out)


if __name__ == "__main__":
    main(
        snakemake.input.fname_snvs_vcf,
        snakemake.params.path_vcf_annotator,
        snakemake.input.fname_genbank_file,
        snakemake.output.fname_snvs_vcf,
    )
