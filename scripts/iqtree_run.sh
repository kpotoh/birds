#!/bin/bash

mkdir -p iqtree_runs/drun1
cd iqtree_runs/drun1
cp ../../trimed_aln_devilworm_clean/* .
cp ../../scheme_devilworm.nex .


# gene partition anc
iqtree2 -p scheme.nex -m MFP+MERGE -asr -nt 5 --prefix label


# KG advise
iqtree2 -s ../folder -m GTR+FO+R6+I -asr -nt AUTO --prefix mephisto_anc
