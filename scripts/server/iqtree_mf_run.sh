#!/bin/bash
#PBS -d .
#PBS -l walltime=50:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=mf

# SCHEME=scheme_birds_max.nex
ALIGNMENT_DIR=aln
CONSTRAINT=constraint.tre

cd ./birds/brun3

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO start building tree" | telegram-send --stdin

$IQTREE -s $ALIGNMENT_DIR -g $CONSTRAINT -m MF -nt $THREADS --prefix $PREFIX
echo "$(date) INFO tree builded" | telegram-send --stdin

