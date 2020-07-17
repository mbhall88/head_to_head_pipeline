import pytest
from concordance import *

from unittest.mock import patch


class TestClassify:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_positionsDontMatch_raisesError(self, mocked_avariant, mocked_bvariant):
        classifier = Classifier()
        mocked_avariant.POS = 1
        mocked_bvariant.POS = 2

        with pytest.raises(IndexError):
            classifier.classify(mocked_avariant, mocked_bvariant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_positionInMask_returnsMasked(self, mocked_variant):
        mask = Bed()
        pos = 2
        mask.positions = set([pos])
        classifier = Classifier(mask=mask)
        mocked_variant.POS = pos

        actual = classifier.classify(mocked_variant, mocked_variant)
        expected = Outcome.Masked

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_positionHasNullButSkipNull_returnsSkippedNull(
        self, mocked_avariant, mocked_bvariant
    ):
        pos = 2
        classifier = Classifier(skip_null=True)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[-1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.SkippedNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothHaveNull_returnsTrueNull(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[-1, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[-1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.TrueNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bHasNullOnly_returnsFalseNull(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[-1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.FalseNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothRef_returnsTrueRef(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, False]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.TrueRef

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bIsRef_returnsFalseRef(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, False]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.FalseRef

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aIsRefBIsAlt_returnsFalseAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[3]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.FalseAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothAlt_returnsTrueAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, 1]]
        mocked_avariant.ALT = ["C"]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[1]]
        mocked_bvariant.ALT = ["C"]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.TrueAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothAltButDifferent_returnsDiffAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, 1]]
        mocked_avariant.ALT = ["C"]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[1]]
        mocked_bvariant.ALT = ["A"]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.DiffAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothFailFilter_returnsBothFailFilter(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = 'b1'
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = 'f0.90;z'

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.BothFailFilter

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aFailFilter_returnsAFailFilter(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = 'b1'
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = None

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.AFailFilter

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bFailFilter_returnsBFailFilter(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = None
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = "foo;bar"

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.BFailFilter

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothHet_returnsBothHet(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, 1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, 1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.BothHet

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aIsHet_returnsAHet(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, 1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.AHet

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bIsHet_returnsBHet(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, 1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.BHet

        assert actual == expected
