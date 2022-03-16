#!/bin/bash
#PBS -d .
#PBS -l walltime=50:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=phylo_ex_constraint

SCHEME=phylo.best_model.nex  # precalculated model on previous run

cd ./birds/brun3

source /home/kpotoh/tools/python_env/bin/activate
echo "$(date) INFO $PREFIX start building tree" | telegram-send --stdin

$IQTREE -p $SCHEME -nt $THREADS --prefix $PREFIX
echo "$(date) INFO $PREFIX tree builded" | telegram-send --stdin