#!/bin/bash

# Mergus_squamatus|Agapornis_pullarius - for birds
# Angiostrongylus_vasorum|Dictyocaulus_eckerti|Oxyuris_equi|Passalurus_ambiguus - for nematoda

LABEL=nematoda

GENCODE=2
if [[ $LABEL -eq "nematoda" ]]; then GENCODE=5; fi
echo "Label $LABEL"
echo -e "Using Gencode $GENCODE\n"

SPECIES=data/interim/species_$LABEL.txt
USED_SP=logs/used_sp_in_aln_$LABEL.log

LOCAL_TMPDIR=/tmp/alignments_${LABEL}
OUTDIR=data/interim/alignments_${LABEL}_clean

rm -fr $LOCAL_TMPDIR $OUTDIR
mkdir -p $LOCAL_TMPDIR $OUTDIR

# Drop pseudogenes and trash
parallel echo Drop pseudogenes from {/} ';' python scripts/qc_alignment.py {} $LOCAL_TMPDIR/{/} --gencode $GENCODE ::: $@  # --dry-run

for name in `cat $SPECIES`
do 
    echo -n -e "$name\t"
    cat $LOCAL_TMPDIR/* | grep -c $name
done > $USED_SP

# Find species with low number of seqs
drop_sp=""
for SP in `egrep -v '1[012]' $USED_SP | cut -f 1`
do
    drop_sp+="$SP|"
done
drop_sp=${drop_sp%"|"}

if [[ $LABEL -eq "birds" ]]; then drop_sp+="|Agapornis_pullarius"; fi
echo -e "Species with low number of genes:\n$drop_sp\n"

# Drop species with low number of genes
for file in $LOCAL_TMPDIR/*
do
    filename=`basename $file`
    echo "Drop species with low number of genes from $filename"
    outfile=$OUTDIR/$(basename $file)
    awk -v pat="$drop_sp" 'BEGIN {RS=">"} $0 !~ pat {printf ">"$o}' $file | sed "s/>>/>/" > $outfile
done

echo -e "\nClean alignments was written to $OUTDIR/"
