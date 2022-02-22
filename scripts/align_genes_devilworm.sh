#!/bin/bash

# USAGE: bash ./align_genes.sh [FILE...]
# FILE must be nucleotide multiple fasta file with only Gene_name in record header
# macse uses 2 threads

MACSE=/home/kpotoh/tools/macse_v2.06.jar
THREADS=12  # 24 / 2

LABEL=devilworm
GENCODE=5
OUTDIR=data/interim/alignments_$LABEL

mkdir -p $OUTDIR
echo created directory $OUTDIR

echo start parallel computing...
echo
# on 24 threads
parallel --jobs $THREADS java -jar $MACSE -prog alignSequences -gc_def $GENCODE -seq {} -out_AA $OUTDIR/{/.}.faa -out_NT $OUTDIR/{/.}.fna ::: $@
