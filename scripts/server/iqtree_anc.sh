#!/bin/bash
#PBS -d .
#PBS -l walltime=50:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=anc_mf

TREE=../phylo.treefile
SCHEME=../scheme_birds_max.nex

cd ./birds/brun3/anc

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO $PREFIX start reconstruction" | telegram-send --stdin

$IQTREE -p $SCHEME -m MFP+MERGE -asr -te $TREE -nt $THREADS --prefix $PREFIX
echo "$(date) INFO $PREFIX reconstructed" | telegram-send --stdin