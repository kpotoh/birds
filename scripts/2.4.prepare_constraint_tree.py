"""
TODO
+ check orders
+ dict[order: [species]]
+ replace order in readed newick  # orner_name --> (...)
+ update odrer tree
"""

from collections import defaultdict
import re
import warnings

import pandas as pd

PATH_TO_SPECIES = "./data/interim/species_birds_in_tree.txt"
PATH_TO_ORDER_TREE = "./data/interim/orders.tre"
PATH_TO_TAXONOMY = "./data/interim/taxonomy-9.1.csv"
PATH_TO_OUT_CONST_TREE = "./data/interim/constraint.tre"

warnings.filterwarnings("ignore")

genus_updater = dict(
    Melophus="Emberiza",
    Megalaima="Psilopogon",
    Uragus="Carpodacus",
    Vestiaria="Drepanis",
    Padda="Lonchura",
)


def read_tree(path: str) -> str:
    with open(path) as fin:
        tree = fin.read().strip()
    return tree


def write_tree(tree: str, path):
    assert isinstance(tree, str)
    with open(path, "w") as fout:
        fout.write(tree)


def read_species(path: str):
    with open(path) as fin:
        species = list(map(str.strip, fin.readlines()))
    return species


def read_taxonomy(path: str):
    taxa = pd.read_csv(path)
    taxa["Species"] = taxa["Species"].str.replace(" ", "_")
    return taxa


def prepare_order_mapping(species, taxa):
    sp2order = dict(taxa[taxa.Species.isin(species)][["Species", "Order"]].values)
    for sp in set(species).difference(taxa.Species):
        genus = sp.split("_")[0]
        genus = genus_updater.get(genus) or genus

        order = taxa[taxa.Genus == genus].Order.values[0]
        assert isinstance(order, str)
        sp2order[sp] = order

    assert len(sp2order) == len(species)

    order2sp = defaultdict(list)
    for sp, order in sp2order.items():
        order2sp[order].append(sp)

    return order2sp


def add_species2tree(order_tree, order2species) -> str:
    orders_from_tree = set(re.findall("(\w+)", order_tree))
    orders_from_data = set(order2species.keys())

    for ordr in orders_from_tree:
        if ordr in orders_from_data:
            species = order2species[ordr]
            sp_str = ",".join(species)
            clade = "({})".format(sp_str) if len(species) > 1 else sp_str
            order_tree = order_tree.replace(ordr, clade)
        else:
            order_tree = order_tree.replace(ordr, "NNNNN")
            print(f"WARNING! Order {ordr} are redutant!")
    return order_tree


def main():
    species = read_species(PATH_TO_SPECIES)
    taxa = read_taxonomy(PATH_TO_TAXONOMY)
    order_tree = read_tree(PATH_TO_ORDER_TREE)

    order2species = prepare_order_mapping(species, taxa)

    for order in order2species:
        if order not in order_tree:
            print(f"WARNING! Order '{order}' not in order tree")

    constraint_tree = add_species2tree(order_tree, order2species)
    write_tree(constraint_tree, PATH_TO_OUT_CONST_TREE)

    return constraint_tree


if __name__ == "__main__":
    main()
