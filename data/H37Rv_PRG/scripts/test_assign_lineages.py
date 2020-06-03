import pytest
from unittest.mock import patch
from io import StringIO
from assign_lineages import (
    PanelVariant,
    RowError,
    Lineage,
    InvalidLineageString,
    load_panel,
    Classifier,
    Genotype,
)


class TestLineageFromStr:
    def test_emptyString_raisesError(self):
        s = ""
        with pytest.raises(InvalidLineageString):
            Lineage.from_str(s)

    def test_withPrefixLineage_ignoresPrefix(self):
        s = "lineageBovis"

        actual = Lineage.from_str(s)
        expected = Lineage(major="Bovis")

        assert actual == expected

    def test_onlyMajor_minorIsNone(self):
        s = "lineage4"

        actual = Lineage.from_str(s).minor

        assert actual is None

    def test_majorAndMinor_allMinorsLumpedInTogether(self):
        s = "lineage4.1.2.3.4.5"

        actual = Lineage.from_str(s)
        expected = Lineage(major="4", minor="1.2.3.4.5")

        assert actual == expected


class TestPanelVariantFromRow:
    def test_emptyRow_raisesError(self):
        row = ""
        with pytest.raises(RowError):
            PanelVariant.from_row(row)

    def test_rowWithNoFields_raisesError(self):
        row = ",,,,"
        with pytest.raises(RowError):
            PanelVariant.from_row(row)

    def test_rowPositionFieldMissing_raisesError(self):
        row = "foo,,,,"
        with pytest.raises(ValueError):
            PanelVariant.from_row(row)

    def test_geneCoordIsntInt_returnsNone(self):
        row = ",".join(
            [
                "lineage1",
                "615938",
                "foo",
                "G/A",
                "368",
                "GAG/GAA",
                "E/E",
                "Rv0524",
                "hemL",
            ]
        )

        actual = PanelVariant.from_row(row)
        expected = PanelVariant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", None)

        assert actual == expected

    def test_geneCoordIdIsEmpty_returnsNone(self):
        row = ",".join(
            [
                "lineage1",
                "615938",
                "",
                "G/A",
                "368",
                "GAG/GAA",
                "E/E",
                "Rv0524",
                "hemL",
            ]
        )

        actual = PanelVariant.from_row(row)
        expected = PanelVariant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", None)

        assert actual == expected

    def test_locusIdIsEmpty_returnsNone(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "", "hemL",]
        )

        actual = PanelVariant.from_row(row)
        expected = PanelVariant(Lineage("1"), 615938, "G", "A", None, "hemL", 1104)

        assert actual == expected

    def test_geneNameIsEmpty_returnsNone(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        actual = PanelVariant.from_row(row)
        expected = PanelVariant(Lineage("1"), 615938, "G", "A", "foo", None, 1104)

        assert actual == expected

    def test_alleleHasNoSeparator_raisesError(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "GA", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        with pytest.raises(ValueError):
            PanelVariant.from_row(row)

    def test_tabAsDelim_returnsExpected(self):
        delim = "\t"
        row = delim.join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        actual = PanelVariant.from_row(row, delim=delim)
        expected = PanelVariant(Lineage("1"), 615938, "G", "A", "foo", None, 1104)

        assert actual == expected


class TestLoadPanel:
    def test_emptyFile_returnsEmpty(self):
        stream = StringIO()

        actual = load_panel(stream)
        expected = dict()

        assert actual == expected

    def test_headerOnly_returnsEmpty(self):
        stream = StringIO("header,line\n")

        actual = load_panel(stream, no_header=False)
        expected = dict()

        assert actual == expected

    def test_noHeader_allPanelVariantsInPanel(self):
        stream = StringIO(
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
            "lineage1.1,4404247,1056,G/A,352,CTG/CTA,L/L,Rv3915,-"
        )

        actual = load_panel(stream)
        expected = {
            615938: PanelVariant(
                Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104
            ),
            4404247: PanelVariant(
                Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056
            ),
        }

        assert actual == expected

    def test_blankLineInFile_noError(self):
        stream = StringIO(
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
            "\n"
            "lineage1.1,4404247,1056,G/A,352,CTG/CTA,L/L,Rv3915,-"
        )

        actual = load_panel(stream)
        expected = {
            615938: PanelVariant(
                Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104
            ),
            4404247: PanelVariant(
                Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056
            ),
        }

        assert actual == expected

    def test_duplicatePosition_raisesError(self):
        stream = StringIO(
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
        )

        with pytest.raises(IndexError):
            load_panel(stream)

    def test_nonDefaultDelim_loadsOk(self):
        stream = StringIO(
            "lineage1|615938|1104|G/A|368|GAG/GAA|E/E|Rv0524|hemL\n"
            "lineage1.1|4404247|1056|G/A|352|CTG/CTA|L/L|Rv3915|-"
        )
        delim = "|"

        actual = load_panel(stream, delim=delim)
        expected = {
            615938: PanelVariant(
                Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104
            ),
            4404247: PanelVariant(
                Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056
            ),
        }

        assert actual == expected

    def test_fileWithNoIssues_loadsOk(self):
        stream = StringIO(
            "lineage,position,gene_coord,allele_change,codon_number,codon_change,amino_acid_change,locus_id,gene_name\n"
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
            "lineage1.1,4404247,1056,G/A,352,CTG/CTA,L/L,Rv3915,-"
        )

        actual = load_panel(stream, no_header=False)
        expected = {
            615938: PanelVariant(
                Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104
            ),
            4404247: PanelVariant(
                Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056
            ),
        }

        assert actual == expected


class TestClassifierIsVariantValid:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_posNotInPanel_returnsFalse(self, mock_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mock_variant.POS = 4
        mock_variant.REF = "G"

        assert not classifier.is_variant_valid(mock_variant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantRefDoesntMatchPanel_returnsFalse(self, mock_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mock_variant.POS = 1
        mock_variant.REF = "d"

        assert not classifier.is_variant_valid(mock_variant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_panelVariantAltNotInVariantAlts_returnsFalse(self, mock_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mock_variant.POS = 1
        mock_variant.REF = "G"
        mock_variant.ALT = ["T", "C"]

        assert not classifier.is_variant_valid(mock_variant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_panelVariantAltInVariantAlts_returnsTrue(self, mock_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mock_variant.POS = 1
        mock_variant.REF = "G"
        mock_variant.ALT = ["T", "A"]

        assert classifier.is_variant_valid(mock_variant)


class TestGenotypeIsNull:
    def test_isNull_returnsTrue(self):
        gt = Genotype(-1, -1, False)

        assert gt.is_null()

    def test_isHom_returnsFalse(self):
        gt = Genotype(-1, 1, False)

        assert not gt.is_null()


class TestGenotypeIsHom:
    def test_isHom_returnsTrue(self):
        gt = Genotype(0, 0, False)

        assert gt.is_hom()

    def test_isHet_returnsFalse(self):
        gt = Genotype(1, 0, False)

        assert not gt.is_hom()

    def test_isNull_returnsFalse(self):
        gt = Genotype(-1, -1, False)

        assert not gt.is_hom()

    def test_oneAlleleIsNull_returnsTrue(self):
        gt = Genotype(-1, 1, False)

        assert gt.is_hom()


class TestGenotypeIsHet:
    def test_isHet_returnsTrue(self):
        gt = Genotype(2, 1, False)

        assert gt.is_het()

    def test_isHom_returnsFalse(self):
        gt = Genotype(0, 0, False)

        assert not gt.is_het()

    def test_oneAlleleIsNull_returnsFalse(self):
        gt = Genotype(-1, 1, False)

        assert not gt.is_het()

    def test_isNull_returnsFalse(self):
        gt = Genotype(-1, -1, False)

        assert not gt.is_het()


class TestGenotypeIsHomRef:
    def test_isHomRef_returnsTrue(self):
        gt = Genotype(0, 0, False)

        assert gt.is_hom_ref()

    def test_isHomRefOneNull_returnsTrue(self):
        gt = Genotype(-1, 0, False)

        assert gt.is_hom_ref()

    def test_isHomAlt_returnsFalse(self):
        gt = Genotype(2, 2, False)

        assert not gt.is_hom_ref()


class TestGenotypeIsHomAlt:
    def test_isHomAlt_returnsTrue(self):
        gt = Genotype(1, 1, False)

        assert gt.is_hom_alt()

    def test_isHomAltOneNull_returnsTrue(self):
        gt = Genotype(-1, 3, False)

        assert gt.is_hom_alt()

    def test_isHomRef_returnsFalse(self):
        gt = Genotype(0, 0, False)

        assert not gt.is_hom_alt()


class TestGenotypeAltIndex:
    def test_isHomAlt_returnsIndex(self):
        gt = Genotype(1, 1, False)

        actual = gt.alt_index()
        expected = 0

        assert actual == expected

    def test_isHomAltWithOneNull_returnsIndex(self):
        gt = Genotype(-1, 5, False)

        actual = gt.alt_index()
        expected = 4

        assert actual == expected

    def test_isHomRef_returnsNone(self):
        gt = Genotype(0, 0, False)

        assert gt.alt_index() is None

    def test_isHet_returnsNone(self):
        gt = Genotype(0, 3, False)

        assert gt.alt_index() is None


class TestClassifierSamplesWithLineageVariant:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noSampleHasLineageVariant_returnsEmpty(self, mocked_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[0, 0, False], [1, 1, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = []

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_oneSampleHasLineageVariant_returnsOneIndex(self, mocked_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[0, 0, False], [2, 2, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = [1]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_allSamplesHaveLineageVariant_returnsAllIndices(self, mocked_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[2, -1, False], [2, 2, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = [0, 1]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_posNotInPanel_returnsEmpty(self, mocked_variant):
        index = {
            4: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[2, -1, False], [2, 2, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = []

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_hetSampleHasLineageVariantIncludeHet_returnsHetIndex(self, mocked_variant):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index, include_het=True)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[-1, -1, False], [2, 0, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = [1]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_hetSampleDoesntHaveLineageVariantIncludeHet_returnsEmpty(
        self, mocked_variant
    ):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index, include_het=True)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A", "N"]
        mocked_variant.genotypes = [[-1, -1, False], [1, 3, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = []

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_hetSampleHasLineageVariantDontIncludeHet_doesntReturnHetIndex(
        self, mocked_variant
    ):
        index = {
            1: PanelVariant(Lineage("1"), 1, "G", "A"),
        }
        classifier = Classifier(index, include_het=False)
        mocked_variant.POS = 1
        mocked_variant.REF = "G"
        mocked_variant.ALT = ["T", "A"]
        mocked_variant.genotypes = [[2, 2, False], [2, 0, False]]

        actual = classifier.samples_with_lineage_variant(mocked_variant)
        expected = [0]

        assert actual == expected
