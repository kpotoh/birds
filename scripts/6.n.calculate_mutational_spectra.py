from collections import defaultdict
from operator import index
import os
import sys
from queue import Queue
from typing import Iterable
from datetime import datetime

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import (
    extract_ff_codons, extract_syn_codons, node_parent, 
    possible_sbs, possible_codons
)

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
        self.gcode = gcode
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]
        self.ff_codons = extract_ff_codons(self.codontable)

        tree = tree = PhyloTree(path_to_tree, format=1)
        anc = pd.read_csv(path_to_states, sep="\t", comment='#')
        leaves = pd.read_csv(path_to_leaves, sep="\t")
        aln_sizes = leaves.groupby("Node").apply(len)
        assert aln_sizes.nunique() == 1, "uncomplete state table"

        states = pd.concat([anc, leaves]).sort_values(["Node", "Part", "Site"])
        mutations, edge_mutspec, total_nucl_freqs = self.extract_mutspec_from_tree(states, tree)

        os.makedirs(out_dir)
        print(f"Output directory {out_dir} created", file=sys.stderr)
        path_to_mutations = os.path.join(out_dir, "mutations.csv")
        path_to_nucl_freqs = os.path.join(out_dir, "nucl_freqs.csv")
        path_to_mutspec = os.path.join(out_dir, "mutspec_{}.csv")
        
        mutations.to_csv(path_to_mutations, index=None)
        total_nucl_freqs.to_csv(path_to_nucl_freqs, index=None)
        for label in self.MUT_LABELS:
            fp_mutspec = path_to_mutspec.format(label)
            edge_mutspec[label].to_csv(fp_mutspec, index=None)

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

    def precalc_node2genome(self, states: pd.DataFrame) -> dict:
        node2genome = dict()
        gr = states.groupby("Node")
        for node in states.Node.unique():
            node2genome[node] = gr.get_group(node).State
        return node2genome

    def extract_mutations(
            self, 
            g1: np.ndarray, g2: np.ndarray, 
            name1: str, name2: str, 
            collect_nucl_freqs=True, context=False,
        ):
        """
        Extract alterations of g2 comparing to g1

        params:
        - g1 - reference genome (parent node)
        - g2 - alternative genome (child node)
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
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"
        if collect_nucl_freqs:
            nucl_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
            codon_freqs = {lbl: defaultdict(int) for lbl in self.MUT_LABELS}
        mutations = []
        for i in range(0, n - 2, 3):
            codon1 = g1[i: i + 3]
            codon2 = g2[i: i + 3]
            codon1_str = "".join(codon1)
            codon2_str = "".join(codon2)
            
            if collect_nucl_freqs:
                for j in range(3):
                    nuc1 = codon1[j]
                    up_nuc1 = g1[i + j - 1]
                    down_nuc1 = g1[i + j + 1]
                    context = f"{up_nuc1}{nuc1}{down_nuc1}"

                    nucl_freqs["all"][nuc1] += 1
                    codon_freqs["all"][context] += 1
                    if j == 2 and self.is_four_fold(codon1_str):
                        nucl_freqs["ff"][nuc1] += 1
                        codon_freqs["ff"][context] += 1

                    # TODO count specific nucl_freqs for syn
                    # if (j == 1 or j == 2)??? and ...:
                    #     nucl_freqs["syn"][nuc1] += 1

            # one codon must contain only sbs and it cannot be indel
            if (codon1 == codon2).sum() != 2 or '-' in codon1 or '-' in codon2:
                continue

            label, aa1, aa2 = self.get_mut_label(codon1_str, codon2_str)

            # collect sbs
            for j in range(3):
                nuc1, nuc2 = codon1[j], codon2[j]
                if nuc1 == nuc2:
                    continue
                if label == 1 and j == 2:
                    label = 2 if self.is_four_fold(codon1_str) else label

                up_nuc1 = g1[i + j - 1]
                down_nuc1 = g1[i + j + 1]
                up_nuc2 = g2[i + j - 1]
                down_nuc2 = g2[i + j + 1]
                if up_nuc1 != up_nuc2 or down_nuc1 != down_nuc2 or up_nuc1 == "-" or down_nuc1 == "-":
                    continue

                sbs = {
                    "RefNode": name1,
                    "AltNode": name2,
                    "Mut": f"{nuc1}>{nuc2}",
                    "MutExt": f"{up_nuc1}[{nuc1}>{nuc2}]{down_nuc1}",
                    "Context": f"{up_nuc1}{nuc1}{down_nuc1}",
                    "RefNucl": nuc1,
                    "AltNucl": nuc2,
                    "Label": label,
                    "Pos": i + j + 1,
                    "PosInCodon": j + 1,
                    "RefCodon": codon1_str,
                    "AltCodon": codon2_str,
                    "RefAa": aa1,
                    "AltAa": aa2,
                }
                mutations.append(sbs)
        
        if len(mutations) > n * 0.1:
            print(f"""
            Warning!
            Ref - {name1}
            Alt - {name2}
            Number of mutations between ref and alt genomes are more than 10% ({n * 0.1}) of the genome length - {len(mutations)}""",
                file=sys.stderr
            )
        mut = pd.DataFrame(mutations)
        if collect_nucl_freqs:
            for lbl in self.MUT_LABELS:
                nucl_freqs[lbl] =  {_nucl:  nucl_freqs[lbl][_nucl]  for _nucl  in "ACGT"}
                codon_freqs[lbl] = {_codon: nucl_freqs[lbl][_codon] for _codon in possible_codons}
            return mut, nucl_freqs
        else:
            return mut, None

    def calculate_mutspec(self, mut: pd.DataFrame, nucl_freqs, label: str):
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
        unobserved_sbs = possible_sbs.difference(mutspec.Mut.values)
        for usbs in unobserved_sbs:
            mutspec_appendix.append({"Mut": usbs, "ObsFr": 0})
        mutspec = pd.concat(
            [mutspec, pd.DataFrame(mutspec_appendix)],
            ignore_index=True
        )
        mutspec["RefNuc"] = mutspec.Mut.str.get(0)
        mutspec["AltNuc"] = mutspec.Mut.str.get(2)
        mutspec["Divisor"] = mutspec.RefNuc.map(nucl_freqs)
        mutspec["RawMutSpec"] = mutspec.ObsFr / mutspec.Divisor
        mutspec["MutSpec"] = mutspec["RawMutSpec"] / mutspec["RawMutSpec"].sum()
        return mutspec

    def extract_mutspec_from_tree(self, states, tree):
        nodes = set(states.Node)
        node2genome = self.precalc_node2genome(states)

        discovered_nodes = set()
        discovered_nodes.add(tree.name)
        Q = Queue()
        Q.put(tree)

        edge_mutspec = defaultdict(list)  # all, syn, ff
        mutations = []
        total_nucl_freqs = []  # dict()
        while not Q.empty():
            cur_node = Q.get()
            for child in cur_node.children:
                Q.put(child)

            if cur_node.name not in discovered_nodes:
                discovered_nodes.add(cur_node.name)
                if cur_node.name not in nodes:
                    continue

                # main process starts here
                parent_node = node_parent(cur_node)
                parent_genome = node2genome[parent_node.name]
                child_genome = node2genome[cur_node.name]
                # collect_nucl_freqs = parent_node.name in total_nucl_freqs
                mut, custom_nucl_freqs = self.extract_mutations(
                    parent_genome.values,
                    child_genome.values,
                    parent_node.name,
                    cur_node.name,
                    # collect_nucl_freqs,
                )
                # custom_nucl_freqs = custom_nucl_freqs or total_nucl_freqs[parent_node.name]
                if len(mut) == 0:
                    continue

                mutations.append(mut)

                cur_nucl_freqs = {"node": parent_node.name}
                for i, lbl in enumerate(self.MUT_LABELS):
                    if lbl == "syn":
                        raise NotImplementedError
                    mutspec = self.calculate_mutspec(mut, custom_nucl_freqs[lbl], label=lbl)
                    mutspec["RefNode"] = parent_node.name
                    mutspec["AltNode"] = cur_node.name
                    edge_mutspec[lbl].append(mutspec)
                    for _nucl in "ACGT":
                        cur_nucl_freqs[f"{_nucl}_{lbl}"] = custom_nucl_freqs[lbl][_nucl]

                total_nucl_freqs.append(cur_nucl_freqs)

        mutations = pd.concat(mutations)
        total_nucl_freqs_df = pd.DataFrame(total_nucl_freqs).drop_duplicates()  # TODO rewrite to normal optimal decision
        # edge_mutspec = list(map(pd.concat, edge_mutspec))
        edge_mutspec_df = {lbl: pd.concat(x) for lbl, x in edge_mutspec.items()}
        return mutations, edge_mutspec_df, total_nucl_freqs_df

    @staticmethod
    def get_common_mutspec(edge_mutspec):
        common_mutspec = edge_mutspec.groupby("Mut")[["ObsFr", "RawMutSpec"]].sum()
        common_mutspec["MutSpec"] = common_mutspec["RawMutSpec"] / \
            common_mutspec["RawMutSpec"].sum()
        return common_mutspec


def main():
    path_to_tree =   "./data/interim/iqtree_runs/brun3/anc_kg.treefile"
    path_to_states = "./data/interim/anc_kg_states_birds.tsv"
    path_to_leaves = "./data/interim/leaves_birds_states.tsv"
    out_dir = "./data/processed/birds"
    out_dir = out_dir + "_" + datetime.now().strftime("%d-%m-%y-%H-%M-%S")
    MutSpec(path_to_tree, path_to_states, path_to_leaves, out_dir)


if __name__ == "__main__":
    main()
