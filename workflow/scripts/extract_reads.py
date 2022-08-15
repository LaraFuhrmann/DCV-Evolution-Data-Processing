#! /usr/bin/env python

"""
extract_reads: created by Tim Stuart: https://timoast.github.io/blog/2015-10-12-extractreads/
"""

import pysam


def get_names(names):
    with open(names, "r") as infile:
        n = infile.read().splitlines()
    if "" in n:
        n.remove("")
    return n


def extract_deletion_read_ids(
    fname_bam, deletion_position, deletion_length, chromosome_id
):
    samfile = pysam.AlignmentFile(fname_bam, "rb")
    # check also position before start position for deletion: - 1
    start_0_based = deletion_position - 1
    print('start_0_based',start_0_based)
    end_1_based = deletion_position + deletion_length
    read_ids = []
    for pileupcolumn in samfile.pileup(
        chromosome_id, start_0_based, start_0_based+1, max_depth=100000
    ):
        if pileupcolumn.pos in range(deletion_position-1, end_1_based+1):
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del and not pileupread.is_refskip :
                    print(pileupread.indel)
                    print(pileupread.alignment.mapping_quality,pileupread.alignment.query_name,"Deletion",pileupcolumn.pos)
                    #if pileupread.indel > 1:
                    #print(pileupread.alignment.cigar)
                    read_ids.append(pileupread.alignment.query_name)

    print(set(read_ids))
    print(len(set(read_ids)))

    return set(read_ids)


def extract_reads(read_ids, fname_bam, fname_out):
    bamfile = pysam.AlignmentFile(fname_bam, "rb")
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    header = bamfile.header.copy()
    out = pysam.Samfile(fname_out, "wb", header=header)
    for name in read_ids:
        try:
            name_indexed.find(name)
        except KeyError:
            pass
        else:
            iterator = name_indexed.find(name)
            for x in iterator:
                out.write(x)


def main(options):

    read_ids = extract_deletion_read_ids(
        options.bam,
        options.deletion_position,
        options.deletion_length,
        options.chromosome_id,
    )
    extract_reads(read_ids, options.bam, options.out)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Extract reads by read name from bam file")
    parser.add_argument("-b", "--bam", help="bam file", required=True)
    parser.add_argument(
        "-pos", "--deletion_position", help="position of deletion", required=True, type=int
    )
    parser.add_argument(
        "-len", "--deletion_length", help="length of deletion", required=True, type=int
    )
    parser.add_argument("-c", "--chromosome_id", help="chromosome_id", required=True)
    parser.add_argument(
        "-o", "--out", help="file name for extracted alignments", required=True
    )
    options = parser.parse_args()
    main(options)
