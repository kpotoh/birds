#!/bin/bash


# LABEL=birds
# OUTDIR=data/interim/trimed_aln_${LABEL}_clean
# mkdir -p $OUTDIR

# for file in data/interim/trimed_aln_$LABEL/*
# do 
#     outfile=$OUTDIR/$(basename $file) 
#     awk 'BEGIN {RS=">"} ! /Mergus_squamatus/ {printf ">"$o}' $file | sed "s/>>/>/" > $outfile
# done




LABEL=devilworm
OUTDIR=data/interim/trimed_aln_${LABEL}_clean
mkdir -p $OUTDIR

for file in data/interim/trimed_aln_$LABEL/*
do 
    outfile=$OUTDIR/$(basename $file) 
    awk 'BEGIN {RS=">"} ! /Angiostrongylus_vasorum|Dictyocaulus_eckerti|Oxyuris_equi|Passalurus_ambiguus/ {printf ">"$o}' $file | sed "s/>>/>/" > $outfile
done
