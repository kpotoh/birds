#!/bin/bash

GENCODE=./data/external/genetic_code2.txt

parallel echo {/} ';' python scripts/qc_alignment.py {} {//}_clean/{/} --gencode $GENCODE ::: $@  # --dry-run
