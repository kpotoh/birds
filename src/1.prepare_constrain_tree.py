"""
TODO
+ check orders
+ dict[order: [species]]
+ replace order in readed newick  # orner_name --> (...)
- update odrer tree
"""

from collections import defaultdict
import re

import pandas as pd

PATH_TO_DATA = "./data/raw/final_birds_list_with_no_mistakes.csv"
PATH_TO_CONSTRAIN = "./data/external/constraint/orders.tre"


def read_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0).drop("X", axis=1)
    orders = df[["Species.name", "Taxonomy"]]
    orders["Species.name"] = orders["Species.name"].str.replace(" ", "_")
    orders["Taxonomy"] = orders.Taxonomy.apply(eval)
    assert len(set(orders.Taxonomy.apply(lambda x: x[12]))) == 1  # birds
    assert "Passeriformes" in set(orders.Taxonomy.apply(lambda x: x[14]))
    orders["Order"] = orders.Taxonomy.apply(lambda x: x[14])
    return orders


def df2dict(df: pd.DataFrame) -> dict:
    assert "Species.name" in df.columns and "Order" in df.columns
    dct = defaultdict(list)
    for name, order in df[["Species.name", "Order"]].values:
        dct[order].append(name)
    return dct


def read_constraint_tree(path: str) -> str:
    with open(path) as fin:
        tree = fin.read().strip()
    return tree


def add_species2tree(constraint_tree, order2species) -> str:
    orders_from_tree = re.findall("(\w+)", constraint_tree)
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

    # assert set(order2species.keys()).  # TODO extend tree to our data

    updated_tree = add_species2tree(constraint_tree, order2species)

    return


if __name__ == "__main__":
    main()
