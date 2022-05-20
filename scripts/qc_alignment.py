import os 
import re
from sys import stderr

import click
from Bio import SeqIO

from utils import read_start_stop_codons

CUTOFF = 0.1
FS_CHAR = "!"
# PATH_TO_GENCODE2 = "./data/external/genetic_code2.txt"

DESCR = f"""
Drop genes with frameshifts ('!') from fasta.
If '!' is in the end of gene it will be ignored.
Pass input and output fasta alignments to run script.
"""


def is_frameshifted(seq: str) -> bool:
    cond = FS_CHAR in seq
    if cond:
        print("!!! Frameshift detected !!!", file=stderr)
    return cond
    

def is_contains_stopcodons(seq: str, stopcodons: set) -> bool:
    n = len(seq)
    assert n % 3 == 0, "seq length is not divisible by 3"

    for i in range(0, n - 2, 3):
        codon = seq[i: i + 3]
        if codon in stopcodons:
            print("!!! Stopcodon detected !!!", file=stderr)
            return True
    return False


@click.command("frameshift_dropper", help=DESCR)
@click.argument("alignment", required=True, type=click.Path(exists=True))  # path to fasta file, containing aligned gene sequences
@click.argument("out-fasta", required=True, type=click.Path(exists=False, writable=True))  # path to output fasta
@click.option(
    "-g", "--gencode", required=False, default=2, show_default=True, help="genetic code")
@click.option(
    "-c", "--cutoff", required=False, default=CUTOFF, show_default=True, 
    help="share of genome length at the end of gene that will not be considered while dropping"
)
def main(alignment, out_fasta, gencode, cutoff: float):
    _, stopcodons = read_start_stop_codons(gencode)
    seqs = SeqIO.parse(alignment, "fasta")
    clean_seqs = []
    for rec in seqs:
        # delete non digit-letter characters in header
        rec.id = re.sub("[^\w_]*", "", rec.id)
        rec.name = re.sub("[^\w_]*", "", rec.name)
        rec.description = ""

        seq = str(rec.seq).upper()
        n = len(seq)
        cons_end = int(n - n * cutoff)
        cons_end -= cons_end % 3
        cons_seq = seq[:cons_end]

        if not is_frameshifted(cons_seq) and not is_contains_stopcodons(seq, stopcodons):
            clean_seqs.append(rec)
        else:
            print(rec.name, "\n", file=stderr)

    outdir = os.path.dirname(out_fasta)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)
    SeqIO.write(clean_seqs, out_fasta, "fasta-2line")


if __name__ == "__main__":
    main()
