from typing import List


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
    if digestion == "none":
        yield from non_specific_digest(seq, min_len, max_len)
    elif digestion == "semi":
        yield from semi_specific_digest(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            miscleavages,
            methionine_cleavage,
        )
    else:
        yield from full_digest(
            seq,
            min_len,
            max_len,
            pre,
            not_post,
            post,
            miscleavages,
            methionine_cleavage,
        )


def non_specific_digest(seq, min_len, max_len):
    seq_len = len(seq)
    for i in range(seq_len + 1):
        for j in range(i + min_len, min(seq_len + 1, i + max_len + 1)):
            if j <= seq_len:
                yield seq[i:j]


def semi_specific_digest(
    seq: str,
    min_len: int,
    max_len: int,
    pre: List[str],
    not_post: List[str],
    post: List[str],
    miscleavages: int,
    methionine_cleavage: bool,
):
    seq_len, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"
    length_accepted = lambda x: x >= min_len and x <= max_len

    for i in range(seq_len):
        is_cleavage_site = is_enzymatic(
            seq[min([seq_len - 1, i])],
            seq[min([seq_len - 1, i + 1])],
            pre,
            not_post,
            post,
        )
        is_methionine_cleavage_site = i == 0 and methionine_cleavage
        if i == seq_len - 1 or is_cleavage_site or is_methionine_cleavage_site:
            # peptides with enzymatic C-terminal (both enzymatic and non-enzymatic N-terminal)
            start = starts[0]
            for j in range(start, min([i + 1, seq_len])):
                pep_len = min([i, seq_len - 1]) - j + 1
                if length_accepted(pep_len):
                    yield (seq[j : i + 1])
            starts.append(i + 1)
            methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
            if len(starts) > miscleavages + 1 + methionine_cleaved or i == seq_len:
                starts = starts[1 + methionine_cleaved :]
        else:  # peptides with non enzymatic C-terminal
            for start in starts:
                pep_len = i - start + 1
                if length_accepted(pep_len) and i + 1 not in starts:
                    yield (seq[start : i + 1])


def full_digest(
    seq: str,
    min_len: int,
    max_len: int,
    pre: List[str],
    not_post: List[str],
    post: List[str],
    miscleavages: int,
    methionine_cleavage: bool,
):
    seq_len, starts = len(seq), [0]
    methionine_cleavage = methionine_cleavage and seq[0] == "M"

    check_pre = len(pre) > 0
    check_post = len(post) > 0

    cleavage_sites = [0] if methionine_cleavage else []
    # HACK: inline if statement instead of using is_enzymatic because it is ~20% faster
    cleavage_sites.extend(
        [
            i
            for i in range(seq_len - 1)
            if (
                check_pre
                and seq[i] in pre
                and not seq[min([seq_len - 1, i + 1])] in not_post
            )
            or (check_post and seq[min([seq_len - 1, i + 1])] in post)
        ]
    )
    cleavage_sites.append(seq_len - 1)
    for i in cleavage_sites:
        for start in starts:
            pep_len = i - start + 1
            if min_len <= pep_len <= max_len:
                yield (seq[start : i + 1])
        starts.append(i + 1)
        methionine_cleaved = int(starts[0] == 0 and methionine_cleavage)
        if len(starts) > miscleavages + 1 + methionine_cleaved:
            starts = starts[1 + methionine_cleaved :]


def is_enzymatic(aa1, aa2, pre, not_post, post):
    return (aa1 in pre and aa2 not in not_post) or (aa2 in post)