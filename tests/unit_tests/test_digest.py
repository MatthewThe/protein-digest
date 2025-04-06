from protein_digest import digest


class TestNonSpecificDigest:
    def test_non_specific_digest(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        assert set(digest.non_specific_digest(seq, min_len, max_len)) == set(
            ["ABCDEF", "BCDEFG", "CDEFGH", "ABCDEFG", "BCDEFGH", "ABCDEFGH"]
        )

    def test_non_specific_digest_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        assert set(digest.non_specific_digest(seq, min_len, max_len)) == set(
            ["ABCDEF", "BCDEFG", "CDEFGH", "ABCDEFG", "BCDEFGH"]
        )


class TestSemiSpecificDigest:
    def test_semi_specific_digest_no_cleavage_site(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "ABCDEFGH", "BCDEFGH", "CDEFGH"])

    def test_semi_specific_digest_methionine_cleavage(self):
        seq = "MABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "MABCDE",
                "MABCDEF",
                "MABCDEFG",
                "MABCDEFGH",
                "ABCDEF",
                "ABCDEFG",
                "ABCDEFGH",
                "BCDEFGH",
                "CDEFGH",
            ]
        )

    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_semi_specific_digest_methionine_cleavage_plus_one_miscleavage(self):
        seq = "MABCDEFKKK"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "MABCDE",
                "MABCDEF",
                "MABCDEFK",
                "MABCDEFKK",
                "ABCDEF",
                "ABCDEFK",
                "ABCDEFKK",
                "BCDEFK",
                "BCDEFKK",
                "CDEFKK",
            ]
        )

    def test_semi_specific_digest_no_cleavage_site_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "BCDEFGH", "CDEFGH"])

    def test_semi_specific_digest_no_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "ABCDEFG", "ABCDEFGK", "BCDEFGK", "CDEFGK"])

    def test_semi_specific_digest_one_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "ABCDEF",
                "ABCDEFG",
                "ABCDEFGK",
                "BCDEFGK",
                "CDEFGK",
                "ABCDEFGKX",
                "BCDEFGKX",
                "CDEFGKX",
                "DEFGKX",
            ]
        )

    def test_semi_specific_digest_not_post(self):
        seq = "ABCDEFKPX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.semi_specific_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(
            [
                "ABCDEF",
                "ABCDEFK",
                "ABCDEFKP",
                "ABCDEFKPX",
                "BCDEFKPX",
                "CDEFKPX",
                "DEFKPX",
            ]
        )


class TestFullDigest:
    def test_full_digest_no_cleavage_site(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGH"])

    def test_full_digest_methionine_cleavage(self):
        seq = "MABCDEFGH"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["MABCDEFGH", "ABCDEFGH"])

    # make sure that the methionine cleavage is not counted as a miscleavage
    def test_full_digest_methionine_cleavage_one_miscleavage(self):
        seq = "MABCDEFGHKKK"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["MABCDEFGHK", "MABCDEFGHKK", "ABCDEFGHK", "ABCDEFGHKK"])

    def test_full_digest_no_cleavage_site_max_len(self):
        seq = "ABCDEFGH"
        min_len = 6
        max_len = 7
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set([])

    def test_full_digest_no_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGK"])

    def test_full_digest_one_miscleavage(self):
        seq = "ABCDEFGKX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 1
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFGK", "ABCDEFGKX"])

    def test_full_digest_not_post(self):
        seq = "ABCDEFKPX"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEFKPX"])

    def test_full_digest_post(self):
        seq = "ABCDEFKPXAAA"
        min_len = 6
        max_len = 30
        pre = [""]
        not_post = [""]
        post = ["K"]
        miscleavages = 0
        methionineCleavage = True
        assert set(
            digest.full_digest(
                seq,
                min_len,
                max_len,
                pre,
                not_post,
                post,
                miscleavages,
                methionineCleavage,
            )
        ) == set(["ABCDEF", "KPXAAA"])