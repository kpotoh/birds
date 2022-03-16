#!/bin/bash
#PBS -d .
#PBS -l walltime=50:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=anc_kg

TREE=phylo.treefile
ALN=anc_only

cd ./birds/brun3

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO $PREFIX start reconstruction" | telegram-send --stdin

$IQTREE -s $ALN -m GTR+FO+R6+I -asr -te $TREE -nt $THREADS --prefix $PREFIX
echo "$(date) INFO $PREFIX reconstructed" | telegram-send --stdin