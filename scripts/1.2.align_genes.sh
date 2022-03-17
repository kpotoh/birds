#!/bin/bash

# USAGE: bash ./align_genes.sh [FILE...]
# FILE must be nucleotide multiple fasta file with only Gene_name in record header
# macse uses 2 threads

# modify this for your custom data
LABEL=birds
GENCODE=2

MACSE=/home/kpotoh/tools/macse_v2.06.jar
THREADS=12  # 24 / 2
OUTDIR=data/interim/alignments_$LABEL
mkdir -p $OUTDIR
echo created directory $OUTDIR

echo start parallel computing...
echo
# on 24 threads
parallel --jobs $THREADS java -jar $MACSE -prog alignSequences -gc_def $GENCODE -seq {} -out_AA $OUTDIR/{/.}.faa -out_NT $OUTDIR/{/.}.fna ::: $@  # --dry-run


# on 2 threads
# for file in $@; do
#     echo "Processing $file"
#     java -jar $MACSE -prog alignSequences -gc_def 2 -seq $file -out_AA $file.aa -out_NT $file.nt
#     echo "$file processing complete"
# done
