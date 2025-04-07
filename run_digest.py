from typing import List
from protein_digest import protein_digest


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


peptides = get_digested_peptides(
    seq="ABCDEFGH",
    min_len=6,
    max_len=30,
    digestion="none",
)
print(peptides)

peptides = get_digested_peptides(
    seq="ABCDEFGH",
    min_len=6,
    max_len=30,
    digestion="semi",
)
print(peptides)

peptides = get_digested_peptides(
    seq="MABCDEFGHKKK",
    min_len=6,
    max_len=30,
    digestion="full",
)
print(peptides)
