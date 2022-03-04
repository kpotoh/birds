# Build bird tree

## Environment
- python 3.8+
- MACSE V2.06
- trimAl v1.4.rev22 build[2015-05-21] [docs](http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2)
- iqtree 2.1.3 [docs](http://www.iqtree.org/doc/iqtree-doc.pdf)

Activate python venv
```
python3 -m venv env_birds
source env_birds/bin/activate
pip install -r requirements.txt
```

## Workflow
### 1 Prepare alignment 
1.0 Extract genes for each bird from raw table
```
python3 scripts/extract_genes.py
```

1.1 Extract fasta for each mitochondrial gene from genes table
```
python3 scripts/split_to_genes.py
```

1.2 Align the sequences in genes separately
```
bash scripts/align_genes.sh data/interim/gene_seqs/*  # need to modify to change gencode
```

1.3 Check if all gaps divisible by 3. Run command from file below and look at last column - it must be all zero else there are some bugs in the alignment
```
cat data/interim/alignments_birds/*.fna | egrep -o "\-*" | sort | uniq -c | awk '{print $1 "\t" length($2) "\t" $2 "\t" length($2)%3}' | tee logs/gaps_birds.log
```

1.4 Drop seqs that contains frameshifts and stopcodons in the middle of the genes
```
bash scripts/qc_aln.sh data/interim/alignments_birds/*.fna  # need to modify gencode before run
```

1.5 Trim alignments
```
bash scripts/trim_alignment.sh data/interim/alignments_birds_clean/*.fna
```

1.6 Repeat step 1.3 if needed

1.7 Check occurence of species in alignment
```
for name in `cat data/interim/species.txt`; do echo -n -e "$name\t"; cat data/interim/trimed_aln_birds/* | grep -c $name; done > logs/used_sp_in_aln.log
for name in `cat data/interim/species_nematoda.txt`; do echo -n -e "$name\t"; cat data/interim/trimed_aln_devilworm/* | grep -c $name; done > logs/used_sp_in_aln_devilworm.log
```

1.8 Drop *Mergus_squamatus* because it has only 8 genes after filtration
```
egrep -v '1[012]' logs/used_sp_in_aln_devilworm.log
bash scripts/qc_aln2.sh  # custom file, change in new run
```

1.9 Repeat step 1.7 
```
for name in `cat data/interim/species.txt`; do echo -n -e "$name\t"; cat data/interim/trimed_aln_birds_clean/* | grep -c $name; done | cut -f 2 | sort | uniq
```

1.n ~~Merge alignments of each gene to one fasta and remember positions of the genes~~ (it will be used in 3rd step)
```
TODO
```

### 2 Prepare constraint tree
2.1 **TODO** from [`orders.tre`](data/external/constraint/orders.tre) and species taxonomy [`bla-bla`](data/external/constraint/) **TODO**
```
python3 scripts/prepare_constrain_tree.py
```

### 3 Prepare nexus file TODO
```
cd data/interim/trimed_aln_birds_clean/  # optional
grep -v -m 1 '>' *.fna | awk -F ':' '{ gsub(".fna", "") ; print "charset", $1, "=", $1 ".fna: 1-" length($2) ";"}'
```
And manually add custom features
 
### 4 Run IQTREE
```
# full model search + tree
iqtree2 -p scheme_birds_max.nex -m MFP+MERGE -nt 8 --prefix phylo

# Ancestral state reconstruction TODO
iqtree2 -anc -nt 8 --prefix anc -p ... -s ...

```



## Devilworm name
- *Halicephalobus mephisto*


## Project structure
```
.
├── data
│   ├── external                                    # data from VB paper (TODO add link)
│   │   ├── constraint
│   │   ├── Input_for_iqtree
│   │   └── raw
│   │       ├── best_iqtree_rooted_415specs.newick
│   │       ├── constraint.7z
│   │       ├── fasttree_777specs.newick
│   │       └── Input_for_iqtree.7z
│   ├── interim
│   └── raw                                         # birds reference genomes data from genbank
│       └── final_birds_list_with_no_mistakes.csv.xz
├── nb
│   └── EDA.ipynb
├── Readme.md
└── src
    └── prepare_constrain_tree.py
```

## References
1. [R. Kimball (2019)](https://www.mdpi.com/1424-2818/11/7/109#supplementary) - constraint tree based on orders
2. [V. Burskaia (2021)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8271140/) - tree building procedure
3. [Birds taxa](https://www.worldbirdnames.org/new/ioc-lists/master-list-2/) - most close to Kimball's taxa (v10.1)
4. [MACSE](https://bioweb.supagro.inra.fr/macse/index.php?menu=doc/dochtml) - tool for multiple alignment of coding sequences
5. [Iqtree](http://www.iqtree.org/) - efficient software for phylogenomic inference
6. [Genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=tgencodes#SG1)
7. [trimAl](http://trimal.cgenomics.org/introduction)
