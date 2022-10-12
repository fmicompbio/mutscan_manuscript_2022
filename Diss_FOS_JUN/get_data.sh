#!/bin/bash

mkdir -p FASTQ
# SRR5952429 (trans input 1)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/009/SRR5952429/SRR5952429_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/009/SRR5952429/SRR5952429_2.fastq.gz
# SRR5952430 (trans input 2)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/000/SRR5952430/SRR5952430_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/000/SRR5952430/SRR5952430_2.fastq.gz
# SRR5952431 (trans input 3)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/001/SRR5952431/SRR5952431_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/001/SRR5952431/SRR5952431_2.fastq.gz
# SRR5952432 (trans output 1)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/002/SRR5952432/SRR5952432_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/002/SRR5952432/SRR5952432_2.fastq.gz
# SRR5952433 (trans output 2)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/003/SRR5952433/SRR5952433_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/003/SRR5952433/SRR5952433_2.fastq.gz
# SRR5952434 (trans output 3)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/004/SRR5952434/SRR5952434_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/004/SRR5952434/SRR5952434_2.fastq.gz
