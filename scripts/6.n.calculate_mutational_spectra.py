import os
import sys
from queue import Queue
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from ete3 import PhyloTree

from utils import extract_ff_codons, node_parent, possible_sbs


class MutSpec:
    def __init__(
            self, 
            path_to_tree, 
            path_to_states, 
            path_to_leaves, 
            out_dir,
            gcode=2,
        ):
        self.gcode = gcode
        self.codontable = CodonTable.unambiguous_dna_by_id[gcode]
        self.ff_codons = extract_ff_codons(self.codontable)

        tree = tree = PhyloTree(path_to_tree, format=1)
        anc = pd.read_csv(path_to_states, sep="\t", comment='#')
        leaves = pd.read_csv(path_to_leaves, sep="\t")
        aln_sizes = leaves.groupby("Node").apply(len)
        assert aln_sizes.nunique() == 1, "uncomplete state table"

        states = pd.concat([anc, leaves]).sort_values(["Node", "Part", "Site"])
        mutations, _ = self.extract_mutspec_from_tree(states, tree, "syn")
        
        os.makedirs(out_dir)
        path_to_mutations = os.path.join(out_dir, "mutations.csv")
        mutations.to_csv(path_to_mutations, index=None)

        path_to_mutspec = os.path.join(out_dir, "mutspec_{}.csv")
        for label in ["all", "syn", "ff"]:
            _, edge_mutspec = self.extract_mutspec_from_tree(states, tree, label)
            edge_mutspec.to_csv(path_to_mutspec.format(label), index=None)

    def is_four_fold(self, codon):
        return codon in self.ff_codons

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

    def extract_mutations(self, g1: np.ndarray, g2: np.ndarray, name1: str, name2: str, context=False):
        """
        TODO: work only with changed positions, not all
        table -> indexes of mutated -> extended indexes with codons 

        params:
        - g1 - reference genome (parent node)
        - g2 - alternative genome (child node)
        """
        n, m = len(g1), len(g2)
        assert n == m, f"genomes lengths are not equal: {n} != {m}"
        assert n % 3 == 0, "genomes length must be divisible by 3 (codon structure)"

        mutations = []
        for i in range(0, n - 2, 3):
            codon1 = g1[i: i + 3]
            codon2 = g2[i: i + 3]
            # one codon must contain only sbs and it cannot be indel
            if (codon1 == codon2).sum() != 2 or '-' in codon1 or '-' in codon2:
                continue

            codon1_str = "".join(codon1)
            codon2_str = "".join(codon2)
            label, aa1, aa2 = self.get_mut_label(codon1_str, codon2_str)

            for j in range(3):
                nuc1, nuc2 = codon1[j], codon2[j]
                if nuc1 == nuc2:
                    continue
                if label == 1 and j == 2:
                    label = 2 if self.is_four_fold(codon1_str) else label

                sbs = {
                    "RefNode": name1,
                    "AltNode": name2,
                    "Mut": f"{nuc1}>{nuc2}",
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
        return mut

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

    def extract_mutspec_from_tree(self, states, tree, label: str):
        nodes = set(states.Node)
        node2genome = self.precalc_node2genome(states)

        discovered_nodes = set()
        discovered_nodes.add(tree.name)
        Q = Queue()
        Q.put(tree)

        edge_mutspec = []
        mutations = []
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

                # TODO
                # ERROR must be determined for all, syn, ff SEPARATELY!!!
                nucl_freqs = parent_genome.value_counts().to_dict() 

                mut = self.extract_mutations(
                    parent_genome.values,
                    child_genome.values,
                    parent_node.name,
                    cur_node.name,
                )
                if len(mut) == 0:
                    continue

                mutations.append(mut)

                mutspec = self.calculate_mutspec(mut, nucl_freqs, label=label)
                mutspec["RefNode"] = parent_node.name
                mutspec["AltNode"] = cur_node.name

                edge_mutspec.append(mutspec)
                # break

        mutations = pd.concat(mutations)
        edge_mutspec = pd.concat(edge_mutspec)
        return mutations, edge_mutspec

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
    out_dir = "./data/processed/birds0"
    MutSpec(path_to_tree, path_to_states, path_to_leaves, out_dir)


if __name__ == "__main__":
    main()
