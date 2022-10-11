#!/usr/bin/env python3
"""
Annotate mutation calls using https://github.com/rpetit3/vcf-annotator.
Position of mutation calls are with respect to the modified parental_stock_consensus
and therefore must be mapped back to EB reference space.

"""
from fuc import pyvcf
import subprocess

def transform_to_EB_space(in_vcf, out_vcf):
    """
    adapt positions of mutation calls such that all is wrt the EB-reference space
    parental_stock.conensus.fasta.
    We are adjusting the parental_stock.consensus.bcftools.fasta sequence such
    that it excludes the UTRs.
    This means we exclude positions [0,24] and [9102, 9263] --
    those positions all have less than 9999 coverage.
    --> to map back we need to add 25 to the POS column.
    """
    vf = pyvcf.VcfFrame.from_file(in_vcf)
    vf.df['POS'] = vf.df['POS']+25
    vf.df['CHROM'] = 'AF014388'
    vf.df = vf.df[vf.df['ALT']!='-']
    #vf.df.loc[vf.df.ALT == "-", "ALT"] = ""
    vf.to_file(out_vcf)

def run_vcf_annotator(in_vcf, path_vcf_annotator, genbank_file, out_vcf):

    subprocess.run(
        [
            "python3",
            str(in_vcf),
            str(path_vcf_annotator),
            str(genbank_file),
            "--output",
            out_vcf,
        ],
        check=True,
    )

def main(fname_snv_in, path_vcf_annotator, fname_genbank_file, fname_snv_out):

    fname_snv_temp = str(fname_snv_in).split('.vcf')[0]+'.temp.vcf'
    transform_to_EB_space(fname_snv_in, fname_snv_temp)
    run_vcf_annotator(fname_snv_temp, path_vcf_annotator, fname_genbank_file, fname_snv_out)


if __name__ == "__main__":
    main(
        snakemake.input.fname_snvs_vcf,
        snakemake.input.fname_genbank_file,
        snakemake.params.path_vcf_annotator,
        snakemake.output.fname_snvs_vcf,
    )
