In this directory we check how many reads align to the host.

Use as reference: https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001215.4/

I did it by hand now:

activate environment with bwa:
`conda activate /cluster/work/bewi/members/lfuhrmann/viloca_applications/hiv_clinical/.snakemake/conda/d84f0052e05dd5536d2e07249e7af6d4_`

`bwa index sequence.fasta`

bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

samples/KO/KOc1/raw_data/:
`bwa mem ../../../../sequence.fasta SRR25258500_1.fastq SRR25258500_2.fastq > aln.sam`

`samtools view -F 0x904 -c aln.sam`
49869 aligned reads

echo $(cat SRR25258500_1.fastq|wc -l)/4|bc
4408143

echo $(cat SRR25258500_2.fastq|wc -l)/4|bc
4408143

WT/WTd10/raw_data/

`bwa mem ../../../../sequence.fasta SRR25258490_1.fastq SRR25258490_2.fastq > aln.sam`
samtools view -F 0x904 -c aln.sam
2125

echo $(cat SRR25258490_1.fastq|wc -l)/4|bc
5185145
echo $(cat SRR25258490_2.fastq|wc -l)/4|bc
5185145

OE/OEa10/raw_data/

bwa mem ../../../../sequence.fasta SRR25258487_1.fastq SRR25258487_2.fastq > aln.sam

samtools view -F 0x904 -c aln.sam
63589

echo $(cat SRR25258487_1.fastq|wc -l)/4|bc
4634873

echo $(cat SRR25258487_2.fastq|wc -l)/4|bc
4634873

## --------
conda activate /cluster/work/bewi/members/lfuhrmann/viloca_applications/hiv_clinical/.snakemake/conda/d84f0052e05dd5536d2e07249e7af6d4_
bwa index GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna

cd samples/KO/KOc1/raw_data/
bwa mem ../../../../GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna SRR25258500_1.fastq SRR25258500_2.fastq > aln.all.sam
samtools view -F 0x904 -c aln.all.sam
396814

cd WT/WTd10/raw_data/
bwa mem ../../../../GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna SRR25258490_1.fastq SRR25258490_2.fastq > aln.all.sam
samtools view -F 0x904 -c aln.all.sam
456489

cd OE/OEa10/raw_data/
bwa mem ../../../../GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna SRR25258487_1.fastq SRR25258487_2.fastq > aln.sam
samtools view -F 0x904 -c aln.sam
614537
