from Bio import SeqIO
from collections import namedtuple
from dataclasses import dataclass
import os
import gzip
import sys
import pandas as pd

### --------------------------------------
### ----- From shorah_snv.py -------------
### --------------------------------------

SNP_id = namedtuple('SNP_id', ['pos', 'var'])
@dataclass
class SNV:
    chrom: str
    pos: int
    ref: str
    var: str
    freq: float = 0.0
    support: float = 0.0

def deletion_length(seq):
    """Determines the length of the deletion. Note that a sequence migth have
       more than one deletion
       seq: substring of the reconstructed haplotype
    """
    count = 0
    for c in seq:
        if c == '-':
            count += 1
        else:
            break
    return count

def parseWindow(shorah_directory, line, ref1, threshold=0.9):
    """SNVs from individual support files, getSNV will build
        the consensus SNVs
        It returns a dictionary called snp with the following structure
        key:   pos.allele (position on the reference file and mutated base)
        value: reference name, position, reference_base, mutated base,
               average number of reads, posterior times average n of reads
    """
    from Bio import SeqIO
    from re import search

    snp = {}
    reads = 0.0
    winFile, chrom, beg, end, cov = line.rstrip().split('\t')
    del([winFile, cov])
    filename = 'w-%s-%s-%s.reads-support.fas' % (chrom, beg, end)

    # take cares of locations/format of support file
    if os.path.exists(filename):
        pass
    elif os.path.exists(shorah_directory+'support/' + filename):
        filename = shorah_directory + 'support/' + filename
    elif os.path.exists(shorah_directory+'support/' + filename + '.gz'):
        filename = shorah_directory+'support/' + filename + '.gz'
    elif os.path.exists(shorah_directory+ filename + '.gz'):
        filename = shorah_directory+ filename + '.gz'

    try:
        if filename.endswith('.gz'):
            window = gzip.open(
                filename, 'rb' if sys.version_info < (3, 0) else 'rt')
        else:
            window = open(filename, 'r')
    except IOError:
        return snp

    beg = int(beg)
    end = int(end)
    refSlice = ref1[chrom][beg - 1:end]
    max_snv = -1
    # sequences in support file exceeding the posterior threshold
    for s in SeqIO.parse(window, 'fasta'):
        seq = str(s.seq).upper()
        match_obj = search('posterior=(.*)\s*ave_reads=(.*)', s.description)
        post, av = float(match_obj.group(1)), float(match_obj.group(2))
        if post >= threshold:
            reads += av
            pos = beg
            tot_snv = 0
            aux_del = -1
            for idx, v in enumerate(refSlice):  # iterate on the reference
                if v != seq[idx]:  # SNV detected, save it
                    if seq[idx] == '-':
                        # Avoid counting multiple times a long deletion in the
                        # same haplotype
                        if idx > aux_del:
                            tot_snv += 1
                            # Check for gap characters and get the deletion
                            # length
                            del_len = deletion_length(seq[idx:])
                            aux_del = idx + del_len
                            snp_id = SNP_id(pos=pos, var=seq[idx:aux_del])

                            if snp_id in snp:
                                # Aggregate counts for long deletions which
                                # are observed in multiple haplotypes
                                snp[snp_id].freq += av
                                snp[snp_id].support += post * av
                            else:
                                # Comply with the convention to report deletion
                                # in VCF format. Position correspond to the
                                # preceding position w.r.t. the reference
                                # without a deletion
                                pos_prev = pos - 1
                                reference_seq = ref1[chrom][
                                    (pos_prev - 1):(pos_prev + del_len)]
                                snp[snp_id] = SNV(
                                    chrom, pos_prev, reference_seq,
                                    reference_seq[0], av, post * av)
                    else:
                        tot_snv += 1
                        snp_id = SNP_id(pos=pos, var=seq[idx])
                        if snp_id in snp:
                            snp[snp_id].freq += av
                            snp[snp_id].support += post * av
                        else:
                            snp[snp_id] = SNV(
                                chrom, pos, v, seq[idx], av, post * av)
                pos += 1
            if tot_snv > max_snv:
                max_snv = tot_snv

    # normalize
    for k, v in snp.items():
        v.support /= v.freq
        v.freq /= reads

    window.close() # TODO

    return snp

def main(fname_reference, shorah_directory):
    """
    shorah_directory: samples/cell_line/passage/variants/SNVs/REGION_1
    """

    ref = dict([[s.id, str(s.seq).upper()] for s in SeqIO.parse(fname_reference, 'fasta')])

    tmp = []

    with open(shorah_directory + 'coverage.txt') as cov_file:
        for line in cov_file:
            snp = parseWindow(shorah_directory, line, ref, threshold=0.9)
            winFile, chrom, beg, end, cov = line.rstrip().split('\t')

            for SNV_id, val in sorted(snp.items()):
                snv_dict = {'winFile': winFile,
                            'chrom': chrom,
                            'start': beg,
                            'end': end,
                            'coverage': cov,
                            'position': val.pos,
                            'ref': val.ref,
                            'var': val.var,
                            'freq': val.freq,
                            'support': val.support}

                tmp.append(snv_dict)

    pd.DataFrame(tmp).to_csv(shorah_directory+'cooccurring_mutations.csv')
    print('Wrote file ', shorah_directory,'cooccurring_mutations.csv')
