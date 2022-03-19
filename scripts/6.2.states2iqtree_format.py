"""
Leaves genomes must be written directly without partitioning - only positions in gene!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 

"""

import os
import re

import click
from Bio import SeqIO
import pandas as pd
import tqdm

# PATH_TO_ALN_DIR = "./data/interim/trimed_aln_devilworm_clean"
# PATH_TO_SCHEME = "./data/interim/iqtree_runs/drun1/anc.best_scheme.nex"  # nexus
# PATH_TO_OUT_STATES = "./data/interim/leaves_states.tsv"


def get_aln_files(path: str):
    assert os.path.isdir(path), "path is not directory"
    raw_files = os.listdir(path)
    files = set([os.path.join(path, x)
                for x in raw_files if x.endswith(".fna")])
    return files


def load_scheme(path: str) -> dict:
    """
    TODO rebuild func, there may be different schemes
    """
    with open(path) as handle:
        raw_file = handle.read()
    charsets = re.findall("charset\s(\w+)\s?=\s?([\w_\.]+)\s?:.+;", raw_file)
    scheme = {i: os.path.basename(fp) "TODO return list..." for i, (gn, fp) in enumerate(charsets, 1)}
    return scheme


def calculate_alignment_length(scheme: dict, aln_path: str) -> int:
    aln_len = 0
    for _, genes in scheme.items():
        for gene_fp in genes:
            filepath = os.path.join(aln_path, gene_fp)
            with open(filepath) as fin:
                fin.readline()
                seq = fin.readline().strip()
                cur_aln_len = len(seq)
                aln_len += cur_aln_len
    return aln_len


def parse_alignment(files: list, scheme: dict, aln_path: str) -> pd.DataFrame:
    data = []
    for part, genes in tqdm.tqdm(scheme.items(), "Parts"):
        previous_part_len = 0
        for gene in genes:
            filepath = os.path.join(aln_path, gene)
            assert filepath in files, f"cannot find file {filepath} from scheme"
            fasta = SeqIO.parse(filepath, "fasta")
            for rec in fasta:
                node = rec.name
                seq = str(rec.seq)
                for site, state in enumerate(seq, 1):
                    pos_data = {
                        "Node": node,
                        "Part": part,
                        "Gene": gene,
                        "Site": site + previous_part_len,
                        "State": state,
                    }
                    for nucl in "ACGT":
                        pos_data[f"p_{nucl}"] = int(nucl == state)

                    data.append(pos_data)

            previous_part_len += len(seq)

    df = pd.DataFrame(data)
    return df


def complete_unk_data(data: pd.DataFrame, aln_len: int) -> pd.DataFrame:
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

    for _un_node in tqdm.tqdm(uncomplete_nodes, "Completing nodes"):
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


@click.command("formatter", help="reformat alignment to states table")
@click.option("--aln", "aln_path", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--scheme", "scheme_path", required=True, type=click.Path(True), help="path to scheme that contain gene splitting info of alignment")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln_path, scheme_path, out):
    aln_files = get_aln_files(aln_path)
    scheme = load_scheme(scheme_path)
    print(scheme)

    aln_data = parse_alignment(aln_files, scheme, aln_path)
    aln_len = calculate_alignment_length(scheme, aln_path)

    data_full = complete_unk_data(aln_data, aln_len)
    assert len(data_full) % aln_len == 0, "something wrong..."

    data_full.to_csv(out, sep="\t", index=None)


if __name__ == "__main__":
    main()
