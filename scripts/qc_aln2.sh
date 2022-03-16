#!/bin/bash


# LABEL=birds
# OUTDIR=data/interim/trimed_aln_${LABEL}_clean
# mkdir -p $OUTDIR

# for file in data/interim/trimed_aln_$LABEL/*
# do 
#     outfile=$OUTDIR/$(basename $file) 
#     awk 'BEGIN {RS=">"} ! /Mergus_squamatus/ {printf ">"$o}' $file | sed "s/>>/>/" > $outfile
# done


# Mergus_squamatus|Agapornis_pullarius - for birds
# Angiostrongylus_vasorum|Dictyocaulus_eckerti|Oxyuris_equi|Passalurus_ambiguus - for nematoda


LABEL=birds
INDIR=data/interim/trimed_aln_$LABEL
OUTDIR=data/interim/trimed_aln_${LABEL}_clean
# INDIR=data/interim/alignments_birds_clean
# OUTDIR=data/interim/alignments_birds_clean_clean
mkdir -p $OUTDIR

for file in $INDIR/*
do 
    outfile=$OUTDIR/$(basename $file) 
    awk 'BEGIN {RS=">"} ! /Mergus_squamatus|Agapornis_pullarius/ {printf ">"$o}' $file | sed "s/>>/>/" > $outfile
done
