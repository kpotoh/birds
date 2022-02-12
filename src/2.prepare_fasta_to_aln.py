"""
TODO
+ sort dataframe
+ replace \s to _ in names
+ drop ND6
+ write 12 files with nucleotides
"""
import os
import pandas as pd

PATH_TO_DATA = "./data/raw/final_birds_list_with_no_mistakes.csv"
PATH_TO_OUT_DIR = "./data/interim/gene_seqs/"


def read_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0).drop("X", axis=1)
    seqs = df[["Species.name", "Gene.name", "Sequence"]]
    seqs["Species.name"] = seqs["Species.name"].str.replace(" ", "_")
    seqs = seqs.sort_values("Species.name")
    seqs["Gene.name"] = seqs["Gene.name"].str.extract("\[(.+)\]")
    seqs = seqs[seqs["Gene.name"] != "ND6"]
    return seqs


def file_write(data: pd.DataFrame, path_to_out: str):
    genes = data["Gene.name"].unique()
    assert os.path.exists(path_to_out), f"{path_to_out} doesn't exist"
    fouts = [open(f"{path_to_out}{gn}.fasta", "w") for gn in genes]

    for i, gn in enumerate(genes):
        cur_seqs = data[data["Gene.name"] == gn]
        cur_fout = fouts[i]
        for sp_name, seq in cur_seqs[["Species.name", "Sequence"]].values:
            cur_fout.write(f">{sp_name}\n")
            cur_fout.write(f"{seq}\n")

    for fout in fouts:
        fout.close()


def main():
    seqs = read_data(PATH_TO_DATA)
    assert all(seqs["Gene.name"] != "ND6"), "col Gene.name contains gene ND6"

    if not os.path.exists(PATH_TO_OUT_DIR):
        os.mkdir(PATH_TO_OUT_DIR)

    file_write(seqs, PATH_TO_OUT_DIR)


if __name__ == "__main__":
    main()
