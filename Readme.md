# Build bird tree

## Environment
- python 3.8+
- iqtree
- [MACSE](https://bioweb.supagro.inra.fr/macse/)
- ...

Activate python venv
```
python3 -m venv env_birds
source env_birds/bin/activate
pip install -r requirements.txt
```

## Workflow
### 1 Prepare alignment 

1.1 Extract fasta for each mitochondrial gene from raw table
```
python3 scripts/prepare_fasta_to_aln.py
```
1.2 Align the sequences in genes separately
```
bash scripts/align_genes.sh data/interim/gene_seqs/*
```
1.3 Merge alignments of each gene to one fasta and remember positions of the genes (it will be used in 3rd step)
```
TODO
```

### 2 Prepare constraint tree
2.1 **TODO** from [`orders.tre`](data/external/constraint/orders.tre) and species taxonomy [`bla-bla`](data/external/constraint/) **TODO**
```
python3 scripts/prepare_constrain_tree.py
```

### 3 Prepare nexus file TODO
 




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