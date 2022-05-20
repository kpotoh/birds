#!/bin/bash
#PBS -d .
#PBS -l walltime=30:00:00,mem=20gb,ncpus=12

IQTREE=/home/kpotoh/tools/iqtree-2.1.3-Linux/bin/iqtree2
THREADS=12
PREFIX=anc_kg
LABEL=nematoda

TREE=anc.treefile
ALN=alignments_nematoda_clean

cd birds_git/data/interim/nrun2

source /home/kpotoh/tools/python_env/bin/activate
telegram-send "$(date) INFO $PREFIX $LABEL start reconstruction"

$IQTREE -s $ALN -m GTR+FO+R6+I -asr -te $TREE -nt $THREADS --prefix $PREFIX
telegram-send "$(date) INFO $PREFIX $LABEL reconstructed"
