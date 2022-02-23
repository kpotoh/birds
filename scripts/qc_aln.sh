#!/bin/bash

GENCODE=./data/external/genetic_code5.txt

parallel echo {/} ';' python scripts/qc_alignment.py {} {//}_clean/{/} --gencode $GENCODE ::: $@  # --dry-run
