import os
import pandas as pd

PATH_TO_DATA = "./data/raw/final_birds_list_with_no_mistakes.csv"
PATH_TO_OUT = "./data/interim/birds_genes.csv"
DROP_ND6 = True

# columns in the raw dataframe
COLUMN_SP_NAME = "Species.name"
COLUMN_GENE_NAME = "Gene.name"
COLUMN_SEQ = "Sequence"


def read_data(path: str, drop_ND6=True) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    seqs = df[[COLUMN_SP_NAME, COLUMN_GENE_NAME, COLUMN_SEQ]]
    seqs[COLUMN_SP_NAME] = seqs[COLUMN_SP_NAME].str.replace(" ", "_")
    seqs = seqs.sort_values(COLUMN_SP_NAME)
    seqs[COLUMN_GENE_NAME] = seqs[COLUMN_GENE_NAME].str.extract("\[(.+)\]")
    if drop_ND6:
        seqs = seqs[seqs[COLUMN_GENE_NAME] != "ND6"]
    seqs.columns = "Species,Gene,Sequence".split(",")
    return seqs


def main():
    seqs = read_data(PATH_TO_DATA, DROP_ND6)
    assert all(seqs["Gene"] != "ND6"), (
        "column Gene contains gene ND6"
    )
    seqs.to_csv(PATH_TO_OUT, index=None)


if __name__ == "__main__":
    main()
