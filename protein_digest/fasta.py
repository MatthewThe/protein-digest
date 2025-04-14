from typing import List, Optional
import collections
import itertools
import logging
from pathlib import Path

from . import digest, digest_rs


logger = logging.getLogger(__name__)


def parse_until_first_space(fasta_id: str) -> str:
    return fasta_id.split(" ")[0]


def get_peptide_to_protein_map(
    fasta_file: str,
    db: str = "concat",
    min_len: int = 6,
    max_len: int = 52,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 2,
    methionine_cleavage: bool = True,
    special_aas: List[str] = ["K", "R"],
    parse_id=parse_until_first_space,
    backend: str="rust",
):
    peptide_to_protein_map = collections.defaultdict(list)
    get_digested_peptides = digest_rs.get_digested_peptides
    if backend == "python":
        get_digested_peptides = digest.get_digested_peptides

    logger.info(f"Parsing fasta file: {Path(fasta_file).name}")
    for protein_idx, (protein, seq) in enumerate(
        read_fasta(fasta_file, db, parse_id, special_aas=special_aas)
    ):
        if protein_idx % 10000 == 0:
            logger.info(f"Digesting protein {protein_idx}")
        seen_peptides = set()
        for peptide in get_digested_peptides(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            digestion,
            miscleavages,
            methionine_cleavage,
        ):
            # peptide = peptide
            if peptide not in seen_peptides:
                seen_peptides.add(peptide)
            peptide_to_protein_map[peptide].append(protein)

    return peptide_to_protein_map


def read_fasta_maxquant(
    file_path: str,
    db: str = "target",
    parse_id=parse_until_first_space,
    special_aas: Optional[List[str]] = None,
    decoy_prefix: str = "REV__",
):
    if special_aas is None:
        special_aas = ["K", "R"]

    if db not in ["target", "decoy", "concat"]:
        raise ValueError("unknown db mode: %s" % db)

    has_special_aas = len(special_aas) > 0
    name, seq = None, []
    with open(file_path, "r") as fp:
        for line in itertools.chain(fp, [">"]):
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seq = "".join(seq)
                    if db in ["target", "concat"]:
                        yield (name, seq)

                    if db in ["decoy", "concat"]:
                        rev_seq = seq[::-1]
                        if has_special_aas:
                            rev_seq = swap_special_aas(rev_seq, special_aas)
                        yield (decoy_prefix + name, rev_seq)

                if len(line) > 1:
                    name, seq = parse_id(line[1:]), []
            else:
                seq.append(line)


read_fasta = read_fasta_maxquant


def swap_special_aas(seq: str, special_aas: List[str]):
    """Swaps the special AAs with its preceding amino acid, as is done in MaxQuant.

    e.g. special_aas = ['R', 'K'] transforms ABCKDEFRK into ABKCDERKF
    """
    seq = list(seq)
    for i in range(1, len(seq)):
        if seq[i] in special_aas:
            swap_positions(seq, i, i - 1)
    seq = "".join(seq)
    return seq


def swap_positions(seq: str, pos1: int, pos2: int):
    seq[pos1], seq[pos2] = seq[pos2], seq[pos1]