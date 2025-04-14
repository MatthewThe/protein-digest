from typing import List

from . import protein_digest


def get_peptide_to_protein_map(
    fasta_file: str,
    db: str = "concat",
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
    special_aas: List[str] = ["K", "R"],
):
    return protein_digest.get_peptide_to_protein_map(
        fasta_file=fasta_file,
        db=db,
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        digestion=digestion,
        miscleavages=miscleavages,
        methionine_cleavage=methionine_cleavage,
        special_aas=special_aas,
    )