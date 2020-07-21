from unittest.mock import patch

import pandas as pd
from concordance import *
from pytest import raises, approx


class TestClassification:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsNull(self, mocked_variant):
        mocked_variant.genotypes = [[-1]]

        actual = Classification.from_variant(mocked_variant)
        expected = Classification.Null

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsHomRef(self, mocked_variant):
        mocked_variant.genotypes = [[0, 0]]

        actual = Classification.from_variant(mocked_variant)
        expected = Classification.Ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsHet(self, mocked_variant):
        mocked_variant.genotypes = [[1, 0]]

        actual = Classification.from_variant(mocked_variant)
        expected = Classification.Het

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsHomAlt(self, mocked_variant):
        mocked_variant.genotypes = [[1, 1]]

        actual = Classification.from_variant(mocked_variant)
        expected = Classification.Alt

        assert actual == expected


class TestClassify:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_positionsDontMatch_raisesError(self, mocked_avariant, mocked_bvariant):
        classifier = Classifier()
        mocked_avariant.POS = 1
        mocked_bvariant.POS = 2

        with raises(IndexError):
            classifier.classify(mocked_avariant, mocked_bvariant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_positionInMask_returnsMasked(self, mocked_avariant, mocked_bvariant):
        mask = Bed()
        pos = 2
        mask.positions = set([pos])
        classifier = Classifier(mask=mask)
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[-1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Null, Classification.Ref, Outcome.Masked

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aHasNull_returnsNull(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[-1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Null, Classification.Ref, Outcome.Null

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothHaveNull_returnsNull(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[-1, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[-1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Null, Classification.Null, Outcome.Null

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bHasNullOnly_returnsFalseNull(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[-1]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Alt, Classification.Null, Outcome.FalseNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothRef_returnsTrueRef(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, -1]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, False]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Ref, Classification.Ref, Outcome.TrueRef

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bIsRef_returnsFalseRef(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, -1]]
        mocked_avariant.ALT = ["C"]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[0, False]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Alt, Classification.Ref, Outcome.FalseRef

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aIsRefBIsAlt_returnsFalseAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[3]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Ref, Classification.Alt, Outcome.FalseAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothAlt_returnsTrueAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, 1]]
        mocked_avariant.ALT = ["C"]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[1]]
        mocked_bvariant.ALT = ["C"]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Alt, Classification.Alt, Outcome.TrueAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothAltButDifferent_returnsDiffAlt(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier()
        mocked_avariant.POS = pos
        mocked_avariant.genotypes = [[1, 1]]
        mocked_avariant.ALT = ["C"]
        mocked_bvariant.POS = pos
        mocked_bvariant.genotypes = [[1]]
        mocked_bvariant.ALT = ["A"]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Alt, Classification.Alt, Outcome.DiffAlt

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothFailFilter_returnsBothFailFilter(
        self, mocked_avariant, mocked_bvariant
    ):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = "b1"
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = "f0.90;z"
        mocked_bvariant.genotypes = [[0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Ref, Classification.Ref, Outcome.BothFailFilter

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_aFailFilter_returnsAFailFilter(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = "b1"
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = None
        mocked_bvariant.genotypes = [[0, 0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Ref, Classification.Ref, Outcome.AFailFilter

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bFailFilter_returnsBFailFilter(self, mocked_avariant, mocked_bvariant):
        pos = 2
        classifier = Classifier(apply_filter=True)
        mocked_avariant.POS = pos
        mocked_avariant.FILTER = None
        mocked_avariant.genotypes = [[0, 0]]
        mocked_bvariant.POS = pos
        mocked_bvariant.FILTER = "foo;bar"
        mocked_bvariant.genotypes = [[0, 0]]

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Classification.Ref, Classification.Ref, Outcome.BFailFilter

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
        expected = Classification.Het, Classification.Het, Outcome.Het

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
        expected = Classification.Het, Classification.Ref, Outcome.Het

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
        expected = Classification.Ref, Classification.Het, Outcome.Het

        assert actual == expected


def create_df(data) -> pd.DataFrame:
    return pd.DataFrame(data, columns=COLUMNS)


class TestClassifyRecall:
    def test_noAltInA_returnsOne(self):
        data = [[1, Classification.Ref, Classification.Ref, Outcome.TrueRef]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 1.0

        assert actual == approx(expected)

    def test_oneAltNotFoundInB_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Ref, Outcome.FalseRef]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_oneAltButMaskedNotFoundInB_returnsOne(self):
        data = [[1, Classification.Alt, Classification.Ref, Outcome.Masked]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 1.0

        assert actual == approx(expected)

    def test_oneAltFoundInB_returnsOne(self):
        data = [[1, Classification.Alt, Classification.Alt, Outcome.TrueAlt]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 1.0

        assert actual == approx(expected)

    def test_oneAltButDiffAltInB_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Alt, Outcome.DiffAlt]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_oneAltBothFailFilter_returnsOne(self):
        data = [[1, Classification.Alt, Classification.Ref, Outcome.BothFailFilter]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 1.0

        assert actual == approx(expected)

    def test_oneAltAFailFilter_returnsOne(self):
        data = [[1, Classification.Alt, Classification.Ref, Outcome.AFailFilter]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 1.0

        assert actual == approx(expected)

    def test_oneAltBFailFilter_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Alt, Outcome.BFailFilter]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_oneAltBMissingPos_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Missing, Outcome.MissingPos]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_oneAltBIsHet_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Het, Outcome.Het]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_oneAltBIsNull_returnsZero(self):
        data = [[1, Classification.Alt, Classification.Null, Outcome.FalseNull]]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.0

        assert actual == approx(expected)

    def test_twoAltBHasOne_returnsHalf(self):
        data = [
            [1, Classification.Alt, Classification.Null, Outcome.FalseNull],
            [2, Classification.Alt, Classification.Alt, Outcome.TrueAlt],
            [3, Classification.Alt, Classification.Alt, Outcome.Masked],
        ]
        df = create_df(data)
        classifier = Classifier()

        actual = classifier.recall(df)
        expected = 0.5

        assert actual == approx(expected)
