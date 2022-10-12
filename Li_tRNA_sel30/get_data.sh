#!/bin/bash

mkdir -p FASTQ
# SRR6813745, SRR6813746, SRR6813747, SRR6813748 (input)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/005/SRR6813745/SRR6813745_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/005/SRR6813745/SRR6813745_2.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/006/SRR6813746/SRR6813746_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/006/SRR6813746/SRR6813746_2.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/007/SRR6813747/SRR6813747_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/007/SRR6813747/SRR6813747_2.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/008/SRR6813748/SRR6813748_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/008/SRR6813748/SRR6813748_2.fastq.gz
# SRR6813755, SRR6813756 (output A)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/005/SRR6813755/SRR6813755_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/005/SRR6813755/SRR6813755_2.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/006/SRR6813756/SRR6813756_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/006/SRR6813756/SRR6813756_2.fastq.gz
# SRR6813757, SRR6813758 (output B)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/007/SRR6813757/SRR6813757_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/007/SRR6813757/SRR6813757_2.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/008/SRR6813758/SRR6813758_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/008/SRR6813758/SRR6813758_2.fastq.gz
# SRR6813759 (output C)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/009/SRR6813759/SRR6813759_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/009/SRR6813759/SRR6813759_2.fastq.gz
# SRR6813760 (output D)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/000/SRR6813760/SRR6813760_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/000/SRR6813760/SRR6813760_2.fastq.gz
# SRR6813761 (output E)
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/001/SRR6813761/SRR6813761_1.fastq.gz
wget -P FASTQ -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR681/001/SRR6813761/SRR6813761_2.fastq.gz
