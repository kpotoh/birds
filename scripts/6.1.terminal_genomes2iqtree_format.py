"""
Leaves genomes must be written directly without partitioning - only positions in gene!!
Only filenames need for this, but if there are preliminary aln concatenation it need additional step

Internal States must be rewritten to similar format 

"""

import os
import re
import sys
import tempfile
from collections import defaultdict
from typing import Dict, Tuple

import click
import pandas as pd
import tqdm
from Bio import SeqIO

# PATH_TO_ALN_DIR = "./data/interim/trimed_aln_devilworm_clean"
# PATH_TO_SCHEME = "./data/interim/iqtree_runs/drun1/anc.best_scheme.nex"  # nexus
# PATH_TO_OUT_STATES = "./data/interim/leaves_states.tsv"
NGENES = 12


def get_aln_files(path: str):
    assert os.path.isdir(path), "path is not directory"
    raw_files = os.listdir(path)
    files = set(
        [os.path.join(path, x) for x in raw_files if x.endswith(".fna")]
    )
    return files


def load_scheme(path: str) -> Dict[str, str]:
    """
    parse files like scheme_birds_genes.nex (just separated genes)

    return dict(charset_lbl: gene_fp)
    """
    with open(path) as handle:
        raw_file = handle.read()
    charsets = re.findall("charset\s(\w+)\s?=\s?([\w_\.]+)\s?:.+;", raw_file)
    scheme = {i: os.path.basename(fp) for i, (gn, fp) in enumerate(charsets, 1)}
    return scheme


def parse_alignment(files: list, scheme: dict, aln_dir, out_dir) -> Tuple[str, int]:
    """
    read fasta files from scheme with alignments and write states to table

    return states table and full alignment length
    """
    _, tmp_fp = tempfile.mkstemp(".tsv", "leaves_states_", dir=out_dir)
    print(f"Temporary state table file - {tmp_fp}", file=sys.stderr)
    tmp_handle = open(tmp_fp, "w")
    columns = "Node Part Site State p_A p_C p_G p_T".split()
    tmp_handle.write("\t".join(columns) + "\n")
    aln_lens = []
    files = set(files)
    history = defaultdict(list)
    for part, gene_fn in tqdm.tqdm(scheme.items(), "Parts"):
        filepath = os.path.join(aln_dir, gene_fn)
        assert filepath in files, f"cannot find file {filepath} from scheme"
        fasta = SeqIO.parse(filepath, "fasta")
        for rec in fasta:
            node = rec.name
            history[node].append(part)
            seq = str(rec.seq)
            for site, state in enumerate(seq, 1):
                pos_data = [node, str(part), str(site), state]
                for nucl in "ACGT":
                    p = int(nucl == state)
                    pos_data.append(str(p))

                tmp_handle.write("\t".join(pos_data) + "\n")
        aln_lens.append(len(seq))

    _pass = False
    for node, parts in history.items():
        if len(parts) == NGENES:
            if _pass:
                continue
            full_parts = parts.copy()
            _pass = False
        else:
            unseen_parts = set(full_parts).difference(parts)
            for unp in unseen_parts:
                for site in range(1, aln_lens[unp - 1] + 1):
                    pos_data = [node, str(part), str(site), "-", "0", "0", "0", "0"]
                    tmp_handle.write("\t".join(pos_data) + "\n")

    tmp_handle.close()
    full_aln_len = sum(aln_lens)
    return tmp_fp, full_aln_len


def _complete_unk_data(data: pd.DataFrame, aln_len: int) -> pd.DataFrame:
    """
    add gaps to genes positions that not used in the alignment

    NOT EVERY GENE SEQUENCE FOR SPECIES PRESENTED IN ALN FILES. 
    WE DROPED SOME GENES WITH LOW QUALITY

    TODO rewrite to more optimal variant without reaing to memory
    """
    aln_sizes = data.groupby("Node").apply(len)
    uncomplete_nodes = aln_sizes[aln_sizes != aln_len].index.values

    some_node = uncomplete_nodes[0]
    while some_node in uncomplete_nodes:
        some_node = data.Node.sample().values[0]
    full_genome_pos = set(
        map(tuple, data[data.Node == some_node][["Part", "Site"]].values)
    )
    for _un_node in tqdm.tqdm(uncomplete_nodes, "Completing nodes"):
        un_genome_pos = set(
            map(tuple, data[data.Node == _un_node][["Part", "Site"]].values)
        )
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
@click.option("--aln", "aln_dir", required=True, type=click.Path(True), help="path to directory with gene alignment files")
@click.option("--scheme", "scheme_path", required=True, type=click.Path(True), help="path to scheme that contain gene splitting info of alignment")
@click.option("--out", required=True, type=click.Path(writable=True), help="path to output states file (tsv)")
def main(aln_dir, scheme_path, out):
    out_dir = os.path.dirname(out)
    aln_files = get_aln_files(aln_dir)
    scheme = load_scheme(scheme_path)
    print(scheme)

    tmp_aln_fp, aln_len = parse_alignment(aln_files, scheme, aln_dir, out_dir)
    data_full = pd.read_csv(tmp_aln_fp, sep="\t")
    # data_full = complete_unk_data(aln_data, aln_len)
    assert len(data_full) % aln_len == 0, "something wrong..."

    data_full.to_csv(out, sep="\t", index=None)
    os.remove(tmp_aln_fp)


if __name__ == "__main__":
    main()
