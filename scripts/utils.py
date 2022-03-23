import re
from collections import defaultdict
from typing import List, Set, Tuple, Union

from Bio.Data import CodonTable
from Bio.Data.CodonTable import NCBICodonTableDNA

PATH_TO_GENCODE5 = "data/external/genetic_code5.txt"


possible_sbs = {
    'A>C', 'A>G', 'A>T',
    'C>A', 'C>G', 'C>T',
    'G>A', 'G>C', 'G>T',
    'T>A', 'T>C', 'T>G'
}

# TODO replace by biopython codontable
######################################
######################################

def read_gencode(path: str) -> List[Tuple[str]]:
    pattern = re.compile("([ACGT]{3})\s([A-Z\*])\s([A-Za-z]{3})\s(i?)")
    gencode = []
    with open(path) as fin:
        for line in fin:
            for code in pattern.findall(line):
                gencode.append(code)
    return gencode


def read_start_stop_codons(path: str) -> Tuple[Set[str]]:
    gencode = read_gencode(path)
    startcodons = set()
    stopcodons = set()
    for code in gencode:
        if code[-1] == "i":
            startcodons.add(code[0])
        if code[1] == "*":
            stopcodons.add(code[0])
    return startcodons, stopcodons


######################################
######################################


def extract_ff_codons(codontable: Union[NCBICodonTableDNA, int]):
    if isinstance(codontable, NCBICodonTableDNA):
        pass
    elif isinstance(codontable, int):
        codontable = CodonTable.unambiguous_dna_by_id[codontable]
    else:
        ValueError("passed codontable is not appropriate")

    aa2codons = defaultdict(set)
    for codon, aa in codontable.forward_table.items():
        aa2codons[aa].add(codon)

    ff_codons = set()
    for aa, codons in aa2codons.items():
        if len(codons) >= 4:
            interim_dct = defaultdict(set)
            for codon in codons:
                interim_dct[codon[:2]].add(codon)

            for nn in interim_dct:
                if len(interim_dct[nn]) == 4:
                    ff_codons = ff_codons.union(interim_dct[nn])
    return ff_codons


def extract_syn_codons(codontable: Union[NCBICodonTableDNA, int]):
    """
    extract codons containing mutation that are synonymous

    return dict[codon: set[PosInCodon]] 
    """

    return 


def node_parent(node):
    try:
        return next(node.iter_ancestors())
    except BaseException:
        return None





if __name__ == "__main__":
    print(read_start_stop_codons(PATH_TO_GENCODE5))
