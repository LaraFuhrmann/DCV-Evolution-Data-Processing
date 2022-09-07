We used `DCV_EBref.fasta` to align the reads of the parental stock. Based on that alignment the consensus sequence of the parental stock `parental_stock.consensus.bcftools.fasta` was created.

`DCV_EBref.fasta` contains the UTRs that where not amplified in our experiments. Hence only <1000 reads align to regions of the UTRs.

We are adjusting the `parental_stock.consensus.bcftools.fasta` sequence such that it excludes the UTRs.
This means we exclude positions [0,24] and [9102, 9263] -- those positions all have less than 2000 coverage.


```
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

with open('parental_stock.consensus.fasta', 'w') as out_f:
    for record in SeqIO.parse('parental_stock.consensus.bcftools.fasta', 'fasta'):

        new_seq = record.seq

        new_seq = new_seq[25:9102]
        new_seq = str(new_seq)

        my_seqs = SeqRecord(Seq(new_seq), id = "parental_stock_consensus")
        SeqIO.write(my_seqs, out_f, "fasta")

```
