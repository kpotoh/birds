#!/bin/bash

cat data/interim/alignments_birds/*.fna | egrep -o "\-*" | sort | uniq -c | awk '{print $1 "\t" length($2) "\t" $2 "\t" length($2)%3}' | tee logs/gaps_birds.log

cat data/interim/alignments_devilworm/*.fna | egrep -o "\-*" | sort | uniq -c | awk '{print $1 "\t" length($2) "\t" $2 "\t" length($2)%3}' | tee logs/gaps_devilworm.log

