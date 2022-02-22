from typing import List, Tuple, Set
import re

PATH_TO_GENCODE5 = "data/external/genetic_code5.txt"


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


if __name__ == "__main__":
    print(read_start_stop_codons(PATH_TO_GENCODE5))
