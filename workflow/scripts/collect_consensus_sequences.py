from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(fnames_consensus, fname_merged_out):

    with open(fname_merged_out, 'w') as out_f:
        my_seqs = []
        for fname_consensus in fnames_consensus:
            experiment = fname_consensus.split("/")[-6]
            patient = fname_consensus.split("/")[-4]
            date = fname_consensus.split("/")[-3]

            for record in SeqIO.parse(fname_consensus, 'fasta'):

                new_seq = str(record.seq)
                new_id = experiment+'-'+patient+'-'+date
                my_seqs.append(SeqRecord(Seq(new_seq), id = new_id))

        SeqIO.write(my_seqs, out_f, "fasta")


if __name__ == "__main__":
    main(
        snakemake.input.fnames_consensus,
        snakemake.output.fname_merged,
    )
