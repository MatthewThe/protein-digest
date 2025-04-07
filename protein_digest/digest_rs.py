from typing import List

from . import protein_digest


def get_digested_peptides(
    seq: str,
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
):
    return protein_digest.get_digested_peptides(
        seq=seq,
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        digestion=digestion,
        miscleavages=miscleavages,
        methionine_cleavage=methionine_cleavage,
    )


def non_specific_digest(
    seq: str,
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    digestion: str = "full",
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
):
    return protein_digest.get_digested_peptides(
        seq=seq,
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        digestion="none",
        miscleavages=miscleavages,
        methionine_cleavage=methionine_cleavage,
    )


def semi_specific_digest(
    seq: str,
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
):
    return protein_digest.get_digested_peptides(
        seq=seq,
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        digestion="semi",
        miscleavages=miscleavages,
        methionine_cleavage=methionine_cleavage,
    )


def full_digest(
    seq: str,
    min_len: int = 6,
    max_len: int = 50,
    pre: List[str] = ["K", "R"],
    not_post: List[str] = ["P"],
    post: List[str] = [],
    miscleavages: int = 0,
    methionine_cleavage: bool = True,
):
    return protein_digest.get_digested_peptides(
        seq=seq,
        min_len=min_len,
        max_len=max_len,
        pre=pre,
        not_post=not_post,
        post=post,
        digestion="full",
        miscleavages=miscleavages,
        methionine_cleavage=methionine_cleavage,
    )
