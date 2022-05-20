#!/bin/bash

# USAGE: bash trim_alignments.sh [FILE...]
# FILE must be fasta with alignments with only Gene_name in record header

# modify this for your custom data
LABEL=nematoda
GT=0.95

THREADS=12
OUTDIR=data/interim/trimed_aln_$LABEL

rm -rf $OUTDIR
mkdir -p $OUTDIR

echo -e "Output directory: $OUTDIR"
echo "Start parallel computing..."
parallel --jobs $THREADS trimal -in {} -gt $GT '|' python scripts/fasta2fasta_2line.py '>' $OUTDIR/{/} ::: $@  # --dry-run
echo "Done"
