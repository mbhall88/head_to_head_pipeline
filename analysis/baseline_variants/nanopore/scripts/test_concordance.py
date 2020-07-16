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
        mocked_avariant.num_unknown = 1
        mocked_bvariant.POS = pos
        mocked_bvariant.num_unknown = 0

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.SkippedNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothHaveNull_returnsTrueNull(
            self, mocked_avariant, mocked_bvariant
    ):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.num_unknown = 1
        mocked_bvariant.POS = pos
        mocked_bvariant.num_unknown = 1

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.TrueNull

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bHasNullOnly_returnsFalseNull(
            self, mocked_avariant, mocked_bvariant
    ):
        pos = 2
        classifier = Classifier(skip_null=False)
        mocked_avariant.POS = pos
        mocked_avariant.num_unknown = 0
        mocked_bvariant.POS = pos
        mocked_bvariant.num_unknown = 1

        actual = classifier.classify(mocked_avariant, mocked_bvariant)
        expected = Outcome.FalseNull

        assert actual == expected
