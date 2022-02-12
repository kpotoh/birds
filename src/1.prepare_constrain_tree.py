"""
TODO
+ check orders
+ dict[order: [species]]
+ replace order in readed newick  # orner_name --> (...)
- update odrer tree
"""

from collections import defaultdict
import re
import warnings

import pandas as pd

PATH_TO_DATA = "./data/raw/final_birds_list_with_no_mistakes.csv"
PATH_TO_CONSTRAIN = "./data/external/constraint/orders.tre"
PATH_TO_OLD_ORDERS = "./data/external/constraint/LIST_OF_ORDERS.csv"

warnings.filterwarnings("ignore")


def get_order_from_taxa(taxa: list) -> str:
    """ from such list:
    - ['Aves', 'Neognathae', 'Psittaciformes', 'Psittacidae', 'Pyrrhura']

    return *formes, i.e. 'Psittaciformes'
    """
    order_pattern = "formes"
    for t in taxa:
        if t.endswith(order_pattern):
            return t


def read_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0).drop("X", axis=1)
    orders = df[["Species.name", "Taxonomy"]].drop_duplicates("Species.name")
    orders["Species.name"] = orders["Species.name"].str.replace(" ", "_")
    orders["Taxonomy"] = orders.Taxonomy.apply(eval)
    assert len(set(orders.Taxonomy.apply(lambda x: x[12]))) == 1  # birds
    assert "Passeriformes" in set(orders.Taxonomy.apply(lambda x: x[14]))
    orders["Order"] = orders.Taxonomy.apply(get_order_from_taxa)
    orders.drop("Taxonomy", axis=1, inplace=True)
    return orders


def df2dict(df: pd.DataFrame) -> dict:
    assert "Species.name" in df.columns and "Order" in df.columns
    df.drop_duplicates("Species.name", inplace=True)
    dct = defaultdict(list)
    for name, order in df[["Species.name", "Order"]].values:
        dct[order].append(name)
    return dct


def read_constraint_tree(path: str) -> str:
    with open(path) as fin:
        tree = fin.read().strip()
    return tree


def read_old_list_of_orders(path) -> dict:
    sp2order = dict()
    with open(path) as fin:
        for line in fin:
            sp, ordr = line.strip().split(",")
            sp2order[sp] = ordr
    return sp2order


def add_species2tree(constraint_tree, order2species) -> str:
    orders_from_tree = set(re.findall("(\w+)", constraint_tree))
    orders_from_data = set(order2species.keys())

    for ordr in orders_from_tree:
        if ordr in orders_from_data:
            species = order2species[ordr]
            sp_str = ",".join(species)
            clade = "({})".format(sp_str)
            constraint_tree = constraint_tree.replace(ordr, clade)
        else:
            constraint_tree = constraint_tree.replace(ordr, "NNNNN")
    return constraint_tree


def main():
    orders = read_data(PATH_TO_DATA)
    order2species = df2dict(orders)
    constraint_tree = read_constraint_tree(PATH_TO_CONSTRAIN)
    
    old_sp2order = read_old_list_of_orders(PATH_TO_OLD_ORDERS)
    for ordr, sps in order2species.items():
        for sp in sps:
            old_order = old_sp2order.get(sp)
            if old_order is not None and not old_order == ordr:
                print(old_order, ordr, sp)

    # assert set(order2species.keys()).  # TODO extend tree to our data

    updated_tree = add_species2tree(constraint_tree, order2species)

    return


if __name__ == "__main__":
    main()
