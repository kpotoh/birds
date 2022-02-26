import os
import re

from Bio import SeqIO
import numpy as np
import pandas as pd
import tqdm

PATH_TO_ALN_DIR = "./data/interim/trimed_aln_devilworm_clean"
PATH_TO_SCHEME = "./data/interim/iqtree_runs/drun1/anc.best_scheme.nex"  # nexus
PATH_TO_OUT_STATES = "./data/interim/leaves_states.tsv"


def get_aln_files(path: str):
    assert os.path.isdir(path), "path is not directory"
    raw_files = os.listdir(PATH_TO_ALN_DIR)
    files = set([os.path.join(path, x)
                for x in raw_files if x.endswith(".fna")])
    return files


def load_scheme(path: str) -> dict:
    """
    TODO rebuild func, it work with charset names, not with files
    """
    with open(path) as handle:
        raw_file = handle.read()
    charsets = re.findall("charset\s(\w+)\s=", raw_file)
    scheme = {i: x.split("_") for i, x in enumerate(charsets, 1)}
    return scheme


def calculate_alignment_length(scheme: dict) -> int:
    aln_len = 0
    for part, genes in scheme.items():
        for gene in genes:
            filepath = os.path.join(PATH_TO_ALN_DIR, f"{gene}.fna")
            with open(filepath) as fin:
                fin.readline()
                seq = fin.readline().strip()
                cur_aln_len = len(seq)
                aln_len += cur_aln_len
    return aln_len


def parse_alignment(files: list, scheme: dict) -> pd.DataFrame:
    data = []
    for part, genes in tqdm.tqdm(scheme.items(), "Parts"):
        previous_part_len = 0
        for gene in genes:
            filepath = os.path.join(PATH_TO_ALN_DIR, f"{gene}.fna")
            assert filepath in files, f"cannot find file {filepath} from scheme"
            fasta = SeqIO.parse(filepath, "fasta")
            for rec in fasta:
                node = rec.name
                seq = str(rec.seq)
                for site, state in enumerate(seq, 1):
                    pos_data = {
                        "Node": node,
                        "Part": part,
                        "Site": site + previous_part_len,
                        "State": state,
                    }
                    for nucl in "ACGT":
                        pos_data[f"p_{nucl}"] = int(nucl == state)

                    data.append(pos_data)

            previous_part_len = len(seq)

    df = pd.DataFrame(data)
    return df


def complete_unk_data(data: pd.DataFrame, aln_len: int):
    """
    add gaps to genes positions that not used in the alignment
    """
    data = data.copy()
    aln_sizes = data.groupby("Node").apply(len)
    uncomplete_nodes = aln_sizes[aln_sizes != aln_len].index.values

    some_node = uncomplete_nodes[0]
    while some_node in uncomplete_nodes:
        some_node = data.Node.sample().values[0]
    full_genome_pos = set(map(
        tuple, data[data.Node == some_node][["Part", "Site"]].values))
    
    for _un_node in uncomplete_nodes:
        un_genome_pos = set(map(
            tuple, data[data.Node == _un_node][["Part", "Site"]].values))
        
        gappy_pos = full_genome_pos.difference(un_genome_pos)

        appendix = []
        for part, site in gappy_pos:
            cur_data = {
                "Node": _un_node,
                "Part": part,
                "Site": site,
                "State": "-",
            }
            for nucl in "ACGT":
                cur_data[f"p_{nucl}"] = 0

            appendix.append(cur_data)
        
        appendix_df = pd.DataFrame(appendix)
        data = pd.concat([data, appendix_df])

    return data




def main():
    aln_files = get_aln_files(PATH_TO_ALN_DIR)
    scheme = load_scheme(PATH_TO_SCHEME)
    aln_data = parse_alignment(aln_files, scheme)
    aln_len = calculate_alignment_length(scheme)

    data_full = complete_unk_data(aln_data, aln_len)

    # aln_data.to_csv(PATH_TO_OUT_STATES, sep="\t", index=None)


if __name__ == "__main__":
    main()




error
data[data.Node == some_node][["Part", "Site"]].apply(lambda x: f"{x.Part}_{x.Site}", axis=1).value_counts()
