"""
TODO
+ sort dataframe
+ replace \s to _ in names
+ drop ND6
+ write 12 files with nucleotides
"""
import os
import sys

import click
import pandas as pd

PATH_TO_DATA = "./data/interim/birds_genes.csv"
PATH_TO_OUT_DIR = "./data/interim/birds_gene_seqs/"

COLUMN_SP_NAME = "Species"
COLUMN_GENE_NAME = "Gene"
COLUMN_SEQ = "Sequence"

STOPCODONS = {"TAA", "TAG", "AGA", "AGG"}


def read_data(path: str) -> pd.DataFrame:
    seqs = pd.read_csv(path)[[COLUMN_SP_NAME, COLUMN_GENE_NAME, COLUMN_SEQ]]
    seqs[COLUMN_SP_NAME] = seqs[COLUMN_SP_NAME].str.replace(" ", "_")
    seqs = seqs.sort_values(COLUMN_SP_NAME)
    return seqs


def drop_stopcodon(seq: str):
    n = len(seq)
    r = n % 3
    if r == 0:
        last_codon = seq[-3:]
        if last_codon in STOPCODONS:
            clean_seq = seq.rstrip(last_codon)
        else:
            clean_seq = seq  # TODO is this best action??? - YES
            print(f"last codon '{last_codon}' is not stopcodon, but remainder = {r}", file=sys.stderr)
    else:
        stump = seq[-r:]
        _is_part_of_stop = False
        for codon in STOPCODONS:
            if codon.startswith(stump):
                _is_part_of_stop = True
                break
        if not _is_part_of_stop:
            print(f"Stopcodon stump '{stump}' is not part of nt stopcodons", file=sys.stderr)

        clean_seq = seq.rstrip(stump)
    return clean_seq


def filter_stopcodons(seqs: pd.DataFrame, inplace=False):
    """
    drop stop codons or stumps of its 
    (mt genome is specific and can use poly-A tail for stop codons)
    """
    if not inplace:
        seqs = seqs.copy()
    seqs[COLUMN_SEQ] = seqs[COLUMN_SEQ].apply(drop_stopcodon)
    return seqs


def file_write(data: pd.DataFrame, path_to_out: str):
    genes = data[COLUMN_GENE_NAME].unique()
    assert os.path.exists(path_to_out), f"{path_to_out} doesn't exist"
    fouts = [open(os.path.join(path_to_out, f"{gn}.fasta"), "w") for gn in genes]

    for i, gn in enumerate(genes):
        cur_seqs = data[data[COLUMN_GENE_NAME] == gn]
        cur_fout = fouts[i]
        for sp_name, seq in cur_seqs[[COLUMN_SP_NAME, COLUMN_SEQ]].values:
            cur_fout.write(f">{sp_name}\n")
            cur_fout.write(f"{seq}\n")

    for fout in fouts:
        fout.close()


@click.command("preparator", help="Write fasta files for each gene")
@click.option("--genes", required=True, 
    help=f"path to csv file, containing gene sequences for each species, must contain 3 columns ({COLUMN_SP_NAME}, {COLUMN_SP_NAME}, {COLUMN_SEQ})")
@click.option("--outdir", required=True, help="path to output directory, will rewrite its content")
def main(genes: pd.DataFrame, outdir: str):
    seqs = read_data(genes)
    filter_stopcodons(seqs, inplace=True)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    file_write(seqs, outdir)


if __name__ == "__main__":
    main()
