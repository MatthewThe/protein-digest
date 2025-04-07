# from protein_digest import digest
from protein_digest import digest_rs as digest


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

    def test_full_digest_egfr(self):
        """
        from pyteomics import parser

        # Perform tryptic digestion with 2 missed cleavages and methionine cleavage
        peptides = list(parser.cleave(
            sequence,
            rule=parser.expasy_rules['trypsin'],  # Use trypsin cleavage rule
            missed_cleavages=2,                   # Allow 2 missed cleavages
            min_length=6,
            max_length=30
        ))

        Not sure which one is correct, but pyteomics cleaves after the first 
        arginine even though it is followed by a proline (MR|P), leading to the
        following discrepancy:

        Extra items in the left set:
        'MRPSGTAGAALLALLAALCPASRALEEKK'
        Extra items in the right set:
        'PSGTAGAALLALLAALCPASRALEEK'
        'PSGTAGAALLALLAALCPASRALEEKK'
        'PSGTAGAALLALLAALCPASR'
        """
        seq = "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA"
        min_len = 6
        max_len = 30
        pre = ["K", "R"]
        not_post = ["P"]
        post = []
        miscleavages = 2
        methionineCleavage = False
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
        ) == set(
            [
                'MRPSGTAGAALLALLAALCPASRALEEKK',  # see comment above
                "VCNGIGIGEFKDSLSINATNIK",
                "CDPSCPNGSCWGAGEENCQKLTK",
                "ELREATSPK",
                "ESDCLVCRKFR",
                "ACGADSYEMEEDGVR",
                "DSLSINATNIK",
                "KFRDEATCK",
                "ITDFGLAKLLGAEEK",
                "NYDLSFLKTIQEVAGYVLIALNTVER",
                "SPSDCCHNQCAAGCTGPR",
                "ELREATSPKANK",
                "MFNNCEVVLGNLEITYVQRNYDLSFLK",
                "DPQRYLVIQGDER",
                "GNMYYENSYALAVLSNYDANK",
                "GMNYLEDR",
                "VLGSGAFGTVYKGLWIPEGEK",
                # "PSGTAGAALLALLAALCPASRALEEKK",  # see comment above
                "EITGFLLIQAWPENRTDLHAFENLEIIRGR",
                "ETEFKKIK",
                "TDLHAFENLEIIR",
                "LLGAEEKEYHAEGGK",
                "HFKNCTSISGDLHILPVAFR",
                "GLWIPEGEK",
                "IPVAIK",
                "NVSRGR",
                "TPLLSSLSATSNNSTVACIDRNGLQSCPIK",
                "NGLQSCPIK",
                "LLGICLTSTVQLITQLMPFGCLLDYVREHK",
                "CKKCEGPCR",
                "EITGFLLIQAWPENR",
                "ELVEPLTPSGEAPNQALLRILKETEFK",
                "DNIGSQYLLNWCVQIAKGMNYLEDRR",
                "ITDFGLAK",
                "GNMYYENSYALAVLSNYDANKTGLKELPMR",
                "GMNYLEDRR",
                "MHLPSPTDSNFYR",
                # "PSGTAGAALLALLAALCPASRALEEK",  # see comment above
                "GENSCKATGQVCHALCSPEGCWGPEPR",
                "FRDEATCK",
                "EATSPK",
                "EYHAEGGKVPIKWMALESILHR",
                "GMNYLEDRRLVHR",
                "DTCPPLMLYNPTTYQMDVNPEGK",
                "SLKEISDGDVIISGNKNLCYANTINWK",
                "GKSPSDCCHNQCAAGCTGPR",
                "LLGICLTSTVQLITQLMPFGCLLDYVR",
                "GSHQISLDNPDYQQDFFPK",
                "TDLHAFENLEIIRGRTK",
                "LLQERELVEPLTPSGEAPNQALLRILK",
                "NLQEILHGAVR",
                "LFGTSGQKTK",
                "LPQPPICTIDVYMIMVKCWMIDADSRPKFR",
                "KLFGTSGQKTK",
                "DPQRYLVIQGDERMHLPSPTDSNFYR",
                "VPIKWMALESILHR",
                "FRELIIEFSK",
                "GSHQISLDNPDYQQDFFPKEAKPNGIFK",
                "NVSRGRECVDK",
                "CWMIDADSRPK",
                "GRTKQHGQFSLAVVSLNITSLGLR",
                "CEGPCRKVCNGIGIGEFK",
                "EHKDNIGSQYLLNWCVQIAKGMNYLEDR",
                "DSLSINATNIKHFK",
                "EAKPNGIFKGSTAENAEYLR",
                "TDLHAFENLEIIRGR",
                "NYVVTDHGSCVRACGADSYEMEEDGVRK",
                "GKSPSDCCHNQCAAGCTGPRESDCLVCR",
                "DEATCK",
                "RPAGSVQNPVYHNQPLNPAPSR",
                "KVCQGTSNKLTQLGTFEDHFLSLQR",
                "RRHIVR",
                "MRPSGTAGAALLALLAALCPASRALEEK",
                "IPSIATGMVGALLLLLVVALGIGLFMRRR",
                "LLGAEEK",
                "LVHRDLAARNVLVK",
                "TKQHGQFSLAVVSLNITSLGLRSLK",
                "DCVSCR",
                "ETEFKK",
                "IKVLGSGAFGTVYK",
                "ALEEKK",
                "ELPMRNLQEILHGAVR",
                "SLKEISDGDVIISGNK",
                "QHGQFSLAVVSLNITSLGLR",
                "DCVSCRNVSRGR",
                "EHKDNIGSQYLLNWCVQIAK",
                "IPVAIKELREATSPK",
                "NVLVKTPQHVKITDFGLAK",
                "VCQGTSNK",
                "KLFGTSGQK",
                "TIQEVAGYVLIALNTVER",
                "TKIISNR",
                "EAKPNGIFK",
                "NLCYANTINWKKLFGTSGQK",
                "EISDGDVIISGNKNLCYANTINWKK",
                "EYHAEGGK",
                "MRPSGTAGAALLALLAALCPASR",
                "VCNGIGIGEFK",
                "GSTAENAEYLRVAPQSSEFIGA",
                "DNIGSQYLLNWCVQIAKGMNYLEDR",
                "FRELIIEFSKMAR",
                "ELVEPLTPSGEAPNQALLR",
                "EITGFLLIQAWPENRTDLHAFENLEIIR",
                "IICAQQCSGRCRGK",
                "GNMYYENSYALAVLSNYDANKTGLK",
                "CNLLEGEPR",
                "VKIPVAIKELR",
                "ILKETEFKK",
                "LTKIICAQQCSGRCR",
                "TGLKELPMR",
                "DNIGSQYLLNWCVQIAK",
                "VCQGTSNKLTQLGTFEDHFLSLQR",
                "IPLENLQIIR",
                "ELIIEFSK",
                "GRECVDK",
                "DIVSSDFLSNMSMDFQNHLGSCQK",
                "TIQEVAGYVLIALNTVERIPLENLQIIR",
                "ATGQVCHALCSPEGCWGPEPR",
                "EATSPKANKEILDEAYVMASVDNPHVCR",
                "NYVVTDHGSCVRACGADSYEMEEDGVR",
                "IKVLGSGAFGTVYKGLWIPEGEK",
                "VAPQSSEFIGA",
                "TKQHGQFSLAVVSLNITSLGLR",
                "CPRNYVVTDHGSCVR",
                "QHGQFSLAVVSLNITSLGLRSLK",
                "RLLQERELVEPLTPSGEAPNQALLR",
                "NVLVKTPQHVK",
                "KVCQGTSNK",
                "TCPAGVMGENNTLVWK",
                "LPQPPICTIDVYMIMVKCWMIDADSRPK",
                "CWMIDADSRPKFR",
                "NYDLSFLK",
                "ESDCLVCR",
                "GPDNCIQCAHYIDGPHCVK",
                "VLGSGAFGTVYKGLWIPEGEKVK",
                "GERLPQPPICTIDVYMIMVK",
                "RHIVRK",
                "WMALESILHR",
                "CRGKSPSDCCHNQCAAGCTGPR",
                "NLCYANTINWK",
                "FSNNPALCNVESIQWR",
                "ATGQVCHALCSPEGCWGPEPRDCVSCR",
                "LTQLGTFEDHFLSLQR",
                "EFVENSECIQCHPECLPQAMNITCTGR",
                "NGLQSCPIKEDSFLQR",
                "KCEGPCRK",
                "IPSIATGMVGALLLLLVVALGIGLFMRR",
                "GLWIPEGEKVK",
                "LTKIICAQQCSGR",
                "NLCYANTINWKK",
                "YSFGATCVKKCPR",
                "VLGSGAFGTVYK",
                "KCEGPCR",
                "KIKVLGSGAFGTVYK",
                "IICAQQCSGRCR",
                "GDSFTHTPPLDPQELDILKTVK",
                "LVHRDLAAR",
                "MFNNCEVVLGNLEITYVQR",
                "ALEEKKVCQGTSNK",
                "GRECVDKCNLLEGEPR",
                "ESDCLVCRK",
                "KVCNGIGIGEFK",
                "GENSCK",
                "NYVVTDHGSCVR",
                "IISNRGENSCK",
                "IICAQQCSGR",
                "CWMIDADSRPKFRELIIEFSK",
                "HIVRKR",
                "MARDPQR",
                "ACGADSYEMEEDGVRK",
                "NLQEILHGAVRFSNNPALCNVESIQWR",
                "TPLLSSLSATSNNSTVACIDR",
                "CPRNYVVTDHGSCVRACGADSYEMEEDGVR",
                "ELIIEFSKMAR",
                "EDSFLQR",
                "VCNGIGIGEFKDSLSINATNIKHFK",
                "NCTSISGDLHILPVAFR",
                "SPSDCCHNQCAAGCTGPRESDCLVCRK",
                "ELIIEFSKMARDPQR",
                "CEGPCRK",
                "EISDGDVIISGNK",
                "ILKETEFK",
                "EATSPKANK",
                "LLQERELVEPLTPSGEAPNQALLR",
                "YSFGATCVKK",
                "GSTAENAEYLR",
                "LFGTSGQK",
                "LLGAEEKEYHAEGGKVPIK",
                "IPSIATGMVGALLLLLVVALGIGLFMR",
                "RLLQER",
                "EISDGDVIISGNKNLCYANTINWK",
                "TKIISNRGENSCK",
                "DLAARNVLVKTPQHVK",
                "ANKEILDEAYVMASVDNPHVCR",
                "ELVEPLTPSGEAPNQALLRILK",
                "YLVIQGDER",
                "CDPSCPNGSCWGAGEENCQK",
                "DEATCKDTCPPLMLYNPTTYQMDVNPEGK",
                "SPSDCCHNQCAAGCTGPRESDCLVCR",
                "TGLKELPMRNLQEILHGAVR",
                "CEGPCR",
                # ="PSGTAGAALLALLAALCPASR",  # see comment above
                "ACGADSYEMEEDGVRKCK",
                "TPQHVKITDFGLAKLLGAEEK",
                "KVCNGIGIGEFKDSLSINATNIK",
                "LPQPPICTIDVYMIMVK",
                "DCVSCRNVSR",
                "EILDEAYVMASVDNPHVCR",
                "DLAARNVLVK",
                "MARDPQRYLVIQGDER",
                "RLVHRDLAAR",
                "TVKEITGFLLIQAWPENR",
                "TLRRLLQER",
                "TPQHVKITDFGLAK",
                "GLWIPEGEKVKIPVAIK",
                "YLVIQGDERMHLPSPTDSNFYR",
                "VKIPVAIK",
                "LFGTSGQKTKIISNR",
                "KCPRNYVVTDHGSCVR",
                "YSFGATCVK",
                "ITDFGLAKLLGAEEKEYHAEGGK",
                "EYHAEGGKVPIK",
                "GDSFTHTPPLDPQELDILK",
                "IPVAIKELR",
                "ECVDKCNLLEGEPR",
                "TPQHVK",
            ]
        )
