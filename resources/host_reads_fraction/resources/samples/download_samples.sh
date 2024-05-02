#!/bin/bash

module load sra-toolkit/2.10.5

# ------ Download samples
fastq-dump -O OE/OEa10/raw_data --split-e  SRR25258487
fastq-dump -O WT/WTd10/raw_data --split-e  SRR25258490
fastq-dump -O KO/KOc1/raw_data --split-e  SRR25258500	
