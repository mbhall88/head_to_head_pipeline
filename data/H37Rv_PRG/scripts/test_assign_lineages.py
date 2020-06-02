import pytest
from io import StringIO
from assign_lineages import Variant, RowError, Lineage, InvalidLineageString, load_panel


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


class TestVariantFromRow:
    def test_emptyRow_raisesError(self):
        row = ""
        with pytest.raises(RowError):
            Variant.from_row(row)

    def test_rowWithNoFields_raisesError(self):
        row = ",,,,"
        with pytest.raises(RowError):
            Variant.from_row(row)

    def test_rowPositionFieldMissing_raisesError(self):
        row = "foo,,,,"
        with pytest.raises(ValueError):
            Variant.from_row(row)

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

        actual = Variant.from_row(row)
        expected = Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", None)

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

        actual = Variant.from_row(row)
        expected = Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", None)

        assert actual == expected

    def test_locusIdIsEmpty_returnsNone(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "", "hemL",]
        )

        actual = Variant.from_row(row)
        expected = Variant(Lineage("1"), 615938, "G", "A", None, "hemL", 1104)

        assert actual == expected

    def test_geneNameIsEmpty_returnsNone(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        actual = Variant.from_row(row)
        expected = Variant(Lineage("1"), 615938, "G", "A", "foo", None, 1104)

        assert actual == expected

    def test_alleleHasNoSeparator_raisesError(self):
        row = ",".join(
            ["lineage1", "615938", "1104", "GA", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        with pytest.raises(ValueError):
            Variant.from_row(row)

    def test_tabAsDelim_returnsExpected(self):
        delim = "\t"
        row = delim.join(
            ["lineage1", "615938", "1104", "G/A", "368", "GAG/GAA", "E/E", "foo", "",]
        )

        actual = Variant.from_row(row, delim=delim)
        expected = Variant(Lineage("1"), 615938, "G", "A", "foo", None, 1104)

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

    def test_noHeader_allVariantsInPanel(self):
        stream = StringIO(
            "lineage1,615938,1104,G/A,368,GAG/GAA,E/E,Rv0524,hemL\n"
            "lineage1.1,4404247,1056,G/A,352,CTG/CTA,L/L,Rv3915,-"
        )

        actual = load_panel(stream)
        expected = {
            615938: Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104),
            4404247: Variant(Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056),
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
            615938: Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104),
            4404247: Variant(Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056),
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
            615938: Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104),
            4404247: Variant(Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056),
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
            615938: Variant(Lineage("1"), 615938, "G", "A", "Rv0524", "hemL", 1104),
            4404247: Variant(Lineage("1", "1"), 4404247, "G", "A", "Rv3915", "-", 1056),
        }

        assert actual == expected
