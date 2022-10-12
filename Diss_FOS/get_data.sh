#!/bin/bash

mkdir -p FASTQ
# SRR5952435 (cis input 1)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/005/SRR5952435/SRR5952435_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/005/SRR5952435/SRR5952435_2.fastq.gz
# SRR5952436 (cis input 2)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/006/SRR5952436/SRR5952436_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/006/SRR5952436/SRR5952436_2.fastq.gz
# SRR5952437 (cis input 3)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/007/SRR5952437/SRR5952437_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/007/SRR5952437/SRR5952437_2.fastq.gz
# SRR5952438 (cis output 1)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/008/SRR5952438/SRR5952438_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/008/SRR5952438/SRR5952438_2.fastq.gz
# SRR5952439 (cis output 2)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/009/SRR5952439/SRR5952439_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/009/SRR5952439/SRR5952439_2.fastq.gz
# SRR5952440 (cis output 3)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/000/SRR5952440/SRR5952440_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR595/000/SRR5952440/SRR5952440_2.fastq.gz
