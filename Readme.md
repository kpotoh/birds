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

## Usage

**TODO**


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
3. [Birds taxa](https://www.worldbirdnames.org/new/ioc-lists/master-list-2/) - like in Kimball's tree
