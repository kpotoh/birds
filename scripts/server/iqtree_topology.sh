#!/bin/bash
#PBS -d .
#PBS -l walltime=50:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=topology

SCHEME=phylo.best_model.nex
TREE=both.tree

cd ./birds/brun3

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO $PREFIX start reconstruction" | telegram-send --stdin

iqtree2 -p $SCHEME -z $TREE -n 0 -zb 1000 -au -nt $THREADS --prefix $PREFIX
echo "$(date) INFO $PREFIX reconstructed" | telegram-send --stdin
