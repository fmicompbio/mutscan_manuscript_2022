#!/bin/bash

mkdir -p FASTQ
# GSM3666027 SRR8713496 (bTDP1-IN-5A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/006/SRR8713496/SRR8713496_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/006/SRR8713496/SRR8713496_2.fastq.gz
# GSM3666029 SRR8713498 (bTDP1-IN-7A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/008/SRR8713498/SRR8713498_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/008/SRR8713498/SRR8713498_2.fastq.gz
# GSM3666030 SRR8713499 (bTDP1-IN-8A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/009/SRR8713499/SRR8713499_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/009/SRR8713499/SRR8713499_2.fastq.gz
## Outputs (2 technical replicates per sample)
# GSM3666031 SRR8713500 (bTDP1-OU-5A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/000/SRR8713500/SRR8713500_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/000/SRR8713500/SRR8713500_2.fastq.gz
# GSM3666032 SRR8713501 (bTDP1-OU-5B)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/001/SRR8713501/SRR8713501_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/001/SRR8713501/SRR8713501_2.fastq.gz
# GSM3666035 SRR8713504	(bTDP1-OU-7A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/004/SRR8713504/SRR8713504_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/004/SRR8713504/SRR8713504_2.fastq.gz
# GSM3666036 SRR8713505 (bTDP1-OU-7B)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/005/SRR8713505/SRR8713505_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/005/SRR8713505/SRR8713505_2.fastq.gz
# GSM3666037 SRR8713506 (bTDP1-OU-8A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/006/SRR8713506/SRR8713506_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/006/SRR8713506/SRR8713506_2.fastq.gz
# GSM3666038 SRR8713507 (bTDP1-OU-8B)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/007/SRR8713507/SRR8713507_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR871/007/SRR8713507/SRR8713507_2.fastq.gz