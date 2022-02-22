#!/bin/bash

# USAGE: bash ./align_genes.sh [FILE...]
# FILE must be nucleotide multiple fasta file with only Gene_name in record header
# macse uses 2 threads

# modify this for your custom data
LABEL=birds
GT=0.95

THREADS=12
OUTDIR=data/interim/trimed_aln_$LABEL

mkdir -p $OUTDIR
echo created directory $OUTDIR
echo
echo start parallel computing...
parallel --jobs $THREADS trimal -in {} -gt $GT '|' python scripts/fasta2fasta.py '>' $OUTDIR/{/.}.fna ::: $@  # --dry-run
echo Done