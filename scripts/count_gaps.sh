seq 100 | parallel printf {} ';' printf - ';' egrep -c "[ATGC]\-{$'{}'}[ATGC]" $@

grep -o --- for --only-matching

head data/interim/alignments_devilworm/ATP6.fna | egrep -o "[ATGC]\-{1,100}[ATGC]|^\-{1,200}[ATGC]|[ATGC]\-{1,100}$" | wc -l