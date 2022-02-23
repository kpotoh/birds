"""
USAGE: cat sample.fasta | python scripts/fasta2fasta.py > sample.fasta-2line
"""

from sys import stdin, stdout
from Bio import SeqIO


def main():
    seqs = SeqIO.parse(stdin, "fasta")
    SeqIO.write(seqs, stdout, "fasta-2line")


if __name__ == "__main__":
    main()
