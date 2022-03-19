#!/bin/bash
#PBS -d .
#PBS -l walltime=40:00:00,mem=25gb,ncpus=16

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=16
PREFIX=anc_mf2

TREE=../phylo.treefile
SCHEME=../scheme_birds_max_anc.nex

cd ./birds/brun3/anc

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO $PREFIX start reconstruction" | telegram-send --stdin

$IQTREE -p $SCHEME -m MFP+MERGE -asr -te $TREE -nt $THREADS --prefix $PREFIX
echo "$(date) INFO $PREFIX reconstructed" | telegram-send --stdin