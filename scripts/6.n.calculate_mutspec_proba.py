"""
TODO:
- optimal states reading
- optimal mutations and spectra writing to files
- syn support
+ minimal sequence is gene (Part), not genome
+ 192 comp
- probability approach
- 
"""

import logging
import logging.config
import os
import sys
from collections import defaultdict
from datetime import datetime
from queue import Queue
from typing import Iterable

import yaml
import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import (extract_ff_codons, extract_syn_codons, node_parent,
                   possible_codons, possible_sbs12, possible_sbs192)

EPS = 1e-5
DEFAULT_PATH_TO_LOGCONF = "./configs/log_settings.yaml"


with open(os.environ.get("LOG_CONFIG", DEFAULT_PATH_TO_LOGCONF), "rb") as file:
    log_config = yaml.safe_load(file)
    logging.config.dictConfig(log_config)
    del log_config

logger = logging.getLogger('MutSpecCalc')


class MutSpec:
    def __init__(
            self, 
            path_to_tree,
            path_to_states,
            path_to_leaves,
            out_dir,
            gcode=2,
        ):
        self.MUT_LABELS = ["all", "ff"]  # TODO add syn
        self.nucl_order = ["A", "C", "G", "T"]
        self.gcode = gcode
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]
        self.ff_codons = extract_ff_codons(self.codontable)

        tree, anc, leaves = self.read_data(path_to_tree, path_to_states, path_to_leaves)
        self.node2genome = self.precalc_node2genome(anc, leaves)
        self.nodes = set(self.node2genome.keys())
        mutations, edge_mutspec12, edge_mutspec192, total_nucl_freqs = self.extract_mutspec_from_tree(tree)

        os.makedirs(out_dir)
        print(f"Output directory {out_dir} created", file=sys.stderr)
        path_to_mutations = os.path.join(out_dir, "mutations.csv")
        path_to_nucl_freqs = os.path.join(out_dir, "nucl_freqs.csv")
        path_to_mutspec = os.path.join(out_dir, "mutspec{}_{}.csv")
        
        mutations.to_csv(path_to_mutations, index=None)
        total_nucl_freqs.to_csv(path_to_nucl_freqs, index=None)
        for label in self.MUT_LABELS:
            fp_mutspec12 = path_to_mutspec.format(12, label)
            fp_mutspec192 = path_to_mutspec.format(192, label)
            edge_mutspec12[label].to_csv(fp_mutspec12, index=None)
            edge_mutspec192[label].to_csv(fp_mutspec192, index=None)

    @staticmethod
    def precalc_node2genome(anc: pd.DataFrame, leaves: pd.DataFrame) -> dict:
        node2genome = defaultdict(dict)
        for states in [anc, leaves]:
            logger.debug("Sorting states")
            states = states.sort_values(["Node", "Part", "Site"])
            logger.debug("States sorted")
            gr = states.groupby(["Node", "Part"])
            for (node, part), gene_pos_ids in gr.groups.items():
                gene_df = states.loc[gene_pos_ids]
                gene_states = gene_df[["p_A", "p_C", "p_G", "p_T"]].values
                node2genome[node][part] = gene_states
        logger.info("node2genome mapping builded")
        return node2genome

    @staticmethod
    def read_data(path_to_tree, path_to_states, path_to_leaves, states_dtype=np.float16):
        tree = PhyloTree(path_to_tree, format=1)
        logger.info(f"tree loaded, number of leaf nodes: {len(tree)}, total number of nodes: {len(tree.get_cached_content())}")
        dtypes = {
            "p_A": states_dtype, 
            "p_C": states_dtype, 
            "p_G": states_dtype, 
            "p_T": states_dtype,
            "Site": np.int32,
            "Part": np.int8,
        }
        usecols = ["Node", "Part", "Site", "p_A", "p_C", "p_G", "p_T"]
        anc = pd.read_csv(path_to_states, sep="\t", comment='#', usecols=usecols, dtype=dtypes)
        logger.info(f"anc states loaded, shape = {anc.shape}")
        leaves = pd.read_csv(path_to_leaves, sep="\t", usecols=usecols, dtype=dtypes)
        logger.info(f"leaves states loaded, shape = {leaves.shape}")
        aln_sizes = leaves.groupby("Node").apply(len)
        assert aln_sizes.nunique() == 1, "uncomplete leaves state table"
        return tree, anc, leaves

    def extract_mutspec_from_tree(self, tree):
        logger.info("Start extracting mutspec from tree")
        discovered_nodes = set()
        discovered_nodes.add(tree.name)
        Q = Queue()
        Q.put(tree)

        edge_mutspec12 = defaultdict(list)  # all, syn, ff
        edge_mutspec192 = defaultdict(list)
        mutations = []
        total_nucl_freqs = []  # dict()
        while not Q.empty():
            cur_node = Q.get()
            for child in cur_node.children:
                Q.put(child)

            if cur_node.name not in discovered_nodes:
                discovered_nodes.add(cur_node.name)
                if cur_node.name not in self.nodes:
                    continue
                parent_node = node_parent(cur_node)

                # main process starts here
                child_gene  = self.node2genome[cur_node.name]
                parent_gene = self.node2genome[parent_node.name]
                # collect_nucl_freqs = parent_node.name in total_nucl_freqs
                sp_mutations = []

                custom_nucl_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
                codon_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}

                for gene in parent_gene:
                    genome_ref = parent_gene[gene]
                    genome_alt = child_gene[gene]
                    mut, gene_nucl_freqs, gene_codon_freqs = self.extract_mutations(
                        genome_ref, genome_alt,
                        parent_node.name, cur_node.name,
                        gene,
                        # collect_nucl_freqs,
                    )
                    for lbl in self.MUT_LABELS:
                        for nucl, freq in gene_nucl_freqs[lbl].items():
                            custom_nucl_freqs[lbl][nucl] += freq
                        for trinucl, freq in gene_codon_freqs[lbl].items():
                            codon_freqs[lbl][trinucl] += freq

                    # custom_nucl_freqs = custom_nucl_freqs or total_nucl_freqs[parent_node.name]
                    if len(mut) > 0:
                        sp_mutations.append(mut)
                
                if len(sp_mutations) == 0:
                    continue

                sp_mutations_df = pd.concat(sp_mutations)
                mutations.append(sp_mutations_df)
                cur_nucl_freqs = {"node": parent_node.name}
                for lbl in self.MUT_LABELS:
                    for _nucl in "ACGT":
                        cur_nucl_freqs[f"{_nucl}_{lbl}"] = custom_nucl_freqs[lbl][_nucl]
                    if lbl == "syn":
                        raise NotImplementedError

                    mutspec12 = self.calculate_mutspec12(sp_mutations_df, custom_nucl_freqs[lbl], label=lbl)
                    mutspec12["RefNode"] = parent_node.name
                    mutspec12["AltNode"] = cur_node.name
                    mutspec192 = self.calculate_mutspec192(sp_mutations_df, codon_freqs[lbl], label=lbl)
                    mutspec192["RefNode"] = parent_node.name
                    mutspec192["AltNode"] = cur_node.name

                    edge_mutspec12[lbl].append(mutspec12)
                    edge_mutspec192[lbl].append(mutspec192)

                total_nucl_freqs.append(cur_nucl_freqs)

        mutations_df = pd.concat(mutations)
        total_nucl_freqs_df = pd.DataFrame(total_nucl_freqs).drop_duplicates()  # TODO rewrite to normal optimal decision
        edge_mutspec12_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec12.items()}
        edge_mutspec192_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec192.items()}
        return mutations_df, edge_mutspec12_df, edge_mutspec192_df, total_nucl_freqs_df

    def extract_mutations(
            self, 
            g1: np.ndarray, g2: np.ndarray, 
            name1: str, name2: str, 
            gene,
            collect_nucl_freqs=True, context=False,
        ):
        """
        Extract alterations of g2 comparing to g1

        params:
        - g1 - reference sequence (parent node)
        - g2 - alternative sequence (child node)
        - name1 - node name of ref
        - name2 - node name of alt
        - collect_nucl_freqs - 
        - context - TODO

        conditions:
        - in one codon could be only sbs
        - in the context of one mutation couldn't be other sbs
        - indels are not sbs and codons and contexts with sbs are not considered

        return:
        - mut - dataframe of mutations
        - nucl_freqs - dict[lbl: dict[{ACGT}: int]] - nucleotide frequencies for all, syn and ff positions
        """
        logger.debug(f"extraction of mutations from {name1} to {name2}")
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        if collect_nucl_freqs:
            nucl_freqs = {lbl: defaultdict(float) for lbl in self.MUT_LABELS}
            codon_freqs = {lbl: defaultdict(float) for lbl in self.MUT_LABELS}

        mutations = []
        for i in range(3, n - 3, 3):  # TODO go out from codons, pentanucleotides are required AND replace self.codon_iterator to 5nucl format sampling!!
            cdn1_states = g1[i: i + 3]
            cdn2_states = g2[i: i + 3]
            if np.any(cdn1_states.sum(1) == 0) or np.any(cdn2_states.sum(1) == 0):
                continue
            
            for cdn1, cdn1_proba in self.codon_iterator(cdn1_states):
                cdn1_str = "".join(cdn1)
                for cdn2, cdn2_proba in self.codon_iterator(cdn2_states):
                    cdn2_str = "".join(cdn2)

                    if collect_nucl_freqs:
                        for j in range(3):
                            nuc1 = cdn1[j]
                            # TODO add proba diversity here for context
                            up_nuc1 = g1[i + j - 1]
                            down_nuc1 = g1[i + j + 1]
                            context = f"{up_nuc1}{nuc1}{down_nuc1}"

                            nucl_freqs["all"][nuc1] += 1
                            codon_freqs["all"][context] += 1
                            if j == 2 and self.is_four_fold(cdn1_str):
                                nucl_freqs["ff"][nuc1] += 1
                                codon_freqs["ff"][context] += 1
                            # TODO count specific nucl_freqs for syn
                            # if (j == 1 or j == 2)??? and ...:
                            #     nucl_freqs["syn"][nuc1] += 1

                    # each codon must contain only sbs and it cannot be indel
                    if sum([cdn1[x] == cdn2[x] for x in range(3)]) != 2:
                        continue
                    # if '-' in cdn1 or '-' in cdn2:
                    #     # redundant step
                    #     continue

                    label, aa1, aa2 = self.get_mut_label(cdn1_str, cdn2_str)

                    # collect sbs
                    for j in range(3):
                        nuc1, nuc2 = cdn1[j], cdn2[j]
                        if nuc1 == nuc2:
                            continue
                        if label == 1 and j == 2:
                            label = 2 if self.is_four_fold(cdn1_str) else label

                        up_nuc1 = g1[i + j - 1]
                        down_nuc1 = g1[i + j + 1]
                        up_nuc2 = g2[i + j - 1]
                        down_nuc2 = g2[i + j + 1]
                        if up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2 or up_nuc1 == "-" or down_nuc1 == "-":
                            continue

                        sbs = {
                            "RefNode": name1,
                            "AltNode": name2,
                            "Gene": gene,
                            "Mut": f"{nuc1}>{nuc2}",
                            "MutExt": f"{up_nuc1}[{nuc1}>{nuc2}]{down_nuc1}",
                            "Context": f"{up_nuc1}{nuc1}{down_nuc1}",
                            "RefNucl": nuc1,
                            "AltNucl": nuc2,
                            "Label": label,
                            "Pos": i + j + 1,
                            "PosInCodon": j + 1,
                            "RefCodon": cdn1_str,
                            "AltCodon": cdn2_str,
                            "RefAa": aa1,
                            "AltAa": aa2,
                        }
                        mutations.append(sbs)
        
        # if len(mutations) > n * 0.1:
        #     print(f"""
        #     Warning!
        #     Ref - {name1}
        #     Alt - {name2}
        #     Number of mutations between ref and alt genomes are more than 10% ({n * 0.1}) of the genome length - {len(mutations)}""",
        #         file=sys.stderr
        #     )
        mut = pd.DataFrame(mutations)
        if collect_nucl_freqs:
            for lbl in self.MUT_LABELS:
                nucl_freqs[lbl] =  {_nucl:  nucl_freqs[lbl][_nucl]  for _nucl  in "ACGT"}
                codon_freqs[lbl] = {_codon: codon_freqs[lbl][_codon] for _codon in possible_codons}
            return mut, nucl_freqs, codon_freqs
        else:
            return mut, None, None

    @staticmethod
    def calculate_mutspec12(mut: pd.DataFrame, nucl_freqs, label: str):
        cols = ["Label", "Mut"]
        for c in cols:
            assert c in mut.columns, f"Column {c} is not in mut df"

        labels = {"syn", "ff", "all"}
        if isinstance(label, str):
            label = label.lower()
            if label not in labels:
                raise ValueError(f"pass the appropriate label: {labels}")
            if label == "syn":
                label = 1
            elif label == "ff":
                label = 2
            elif label == "all":
                label = 0
        else:
            raise ValueError(f"pass the appropriate label: {labels}")

        mutspec = mut[mut.Label >= label].Mut.value_counts().reset_index()
        mutspec.columns = ["Mut", "ObsFr"]

        mutspec_appendix = []
        unobserved_sbs = possible_sbs12.difference(mutspec.Mut.values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
        mutspec = pd.concat(
            [mutspec, pd.DataFrame(mutspec_appendix)],
            ignore_index=True
        )
        mutspec["RefNuc"] = mutspec.Mut.str.get(0)
        mutspec["Divisor"] = mutspec.RefNuc.map(nucl_freqs)
        mutspec["RawMutSpec"] = mutspec.ObsFr / mutspec.Divisor
        mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
        mutspec.drop("RefNuc", axis=1, inplace=True)
        return mutspec

    @staticmethod
    def calculate_mutspec192(mut: pd.DataFrame, codon_freqs, label: str):
        cols = ["Label", "MutExt"]
        for c in cols:
            assert c in mut.columns, f"Column {c} is not in mut df"

        available_labels = {"syn", "ff", "all"}
        if isinstance(label, str):
            label = label.lower()
            if label not in available_labels:
                raise ValueError(f"pass the appropriate label: {available_labels}")
            if label == "syn":
                label = 1
            elif label == "ff":
                label = 2
            elif label == "all":
                label = 0
        else:
            raise ValueError(f"pass the appropriate label: {available_labels}")

        mutspec = mut[mut.Label >= label].MutExt.value_counts().reset_index()
        mutspec.columns = ["Mut", "ObsFr"]

        mutspec_appendix = []
        unobserved_sbs = possible_sbs192.difference(mutspec.Mut.values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
        mutspec = pd.concat(
            [mutspec, pd.DataFrame(mutspec_appendix)],
            ignore_index=True
        )
        mutspec["Context"] = mutspec.Mut.str.get(0) + mutspec.Mut.str.get(2) + mutspec.Mut.str.get(-1)
        mutspec["Divisor"] = mutspec.Context.map(codon_freqs)
        mutspec["RawMutSpec"] = (mutspec.ObsFr / mutspec.Divisor).fillna(0)
        mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
        mutspec.drop("Context", axis=1, inplace=True)
        return mutspec

    def codon_iterator(self, codon_states: np.ndarray, cutoff=0.01):
        assert codon_states.shape == (3, 4)
        for i, p1 in enumerate(codon_states[0]):
            if p1 < cutoff:
                continue
            for j, p2 in enumerate(codon_states[1]):
                if p2 < cutoff:
                    continue
                for k, p3 in enumerate(codon_states[2]):
                    if p3 < cutoff:
                        continue
                    codon_proba = p1 * p2 * p3
                    if codon_proba < cutoff:
                        continue
                    codon = [self.nucl_order[x] for x in [i, j, k]]
                    yield codon, codon_proba

    def is_four_fold(self, codon):
        return codon in self.ff_codons

    def is_syn(self, codon, pos_in_codon):
        raise NotImplementedError
        # return pos_in_codon in self.syn_codons.get(codon)

    def get_mut_label(self, codon1: str, codon2: str):
        """
        returned labels:
        - -1 - error mutation (contains stopcodon)
        -  0 - usual mutation
        -  1 - synonimous mutation
        """
        assert codon1 != codon2, "codons must be different"
        aa1 = self.codontable.forward_table.get(codon1, "*")
        aa2 = self.codontable.forward_table.get(codon2, "*")
        if aa1 == "*" or aa2 == "*":
            label = -1
        elif aa1 == aa2:
            label = 1
        else:
            label = 0

        return label, aa1, aa2

def main():
    path_to_tree =   "./data/interim/iqtree_runs/brun3/anc_kg.treefile"
    path_to_states = "./data/interim/anc_kg_states_birds.tsv"
    path_to_leaves = "./data/interim/leaves_birds_states.tsv"
    out_dir = "./data/processed/birds"
    out_dir = out_dir + "_" + datetime.now().strftime("%d-%m-%y-%H-%M-%S")
    MutSpec(path_to_tree, path_to_states, path_to_leaves, out_dir)


if __name__ == "__main__":
    main()
