from unittest.mock import patch

import pytest
from apply_variants_to_loci import Record


class TestRecordApplyVariant:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_refDoesntMatch_raisesError(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "hello world"
        record = Record(name, comment, seq)
        relative_start = 1
        mocked_variant.POS = 1
        mocked_variant.REF = "hall"

        with pytest.raises(ValueError):
            record.apply_variant(mocked_variant, relative_start)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noAlt_returnsEmpty(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "hello world"
        record = Record(name, comment, seq)
        relative_start = 1
        mocked_variant.POS = 1
        mocked_variant.REF = "hell"
        mocked_variant.ALT = []

        actual = record.apply_variant(mocked_variant, relative_start)
        expected = []

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_oneAlt_returnsOneRecord(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "hello world"
        record = Record(name, comment, seq)
        relative_start = 1
        mocked_variant.POS = 1
        mocked_variant.REF = "hello"
        mocked_variant.ALT = ["hola"]

        with patch("uuid.uuid4", return_value=name) as mocked_uuid4:
            actual = record.apply_variant(mocked_variant, relative_start)
            assert mocked_uuid4.call_count == 1
        expected_comment = f"POS=1|ALT=0|ALT_POS_IN_SEQ=[0,4)"
        expected = [Record(name, expected_comment, "hola world")]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_twoAlts_returnsTwoRecord(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "hello world"
        record = Record(name, comment, seq)
        relative_start = 2
        max_indel_len = None
        mocked_variant.POS = 4
        mocked_variant.REF = "l"
        mocked_variant.ALT = ["ins", "y"]

        with patch("uuid.uuid4", return_value=name) as mocked_uuid4:
            actual = record.apply_variant(
                mocked_variant, relative_start, max_indel_len=max_indel_len
            )
            assert mocked_uuid4.call_count == 2
        expected_comment1 = f"POS=4|ALT=0|ALT_POS_IN_SEQ=[2,5)"
        expected_comment2 = f"POS=4|ALT=1|ALT_POS_IN_SEQ=[2,3)"
        expected = [
            Record(name, expected_comment1, "heinslo world"),
            Record(name, expected_comment2, "heylo world"),
        ]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_twoAltsOneTooLong_returnsOneRecord(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "hello world"
        record = Record(name, comment, seq)
        max_indel_len = 3
        relative_start = 2
        mocked_variant.POS = 4
        mocked_variant.REF = "l"
        mocked_variant.ALT = ["ins", "longlong"]

        with patch("uuid.uuid4", return_value=name) as mocked_uuid4:
            actual = record.apply_variant(
                mocked_variant, relative_start, max_indel_len=max_indel_len
            )
            assert mocked_uuid4.call_count == 1
        expected_comment1 = f"POS=4|ALT=0|ALT_POS_IN_SEQ=[2,5)"
        expected = [
            Record(name, expected_comment1, "heinslo world"),
        ]

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantStartsInPreviousLoci_skips(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "TACGTfoobar"
        record = Record(name, comment, seq)
        relative_start = 10
        mocked_variant.POS = 7
        mocked_variant.REF = "ACGTACGT"
        mocked_variant.ALT = ["abcdefgh"]

        with patch("uuid.uuid4", return_value=name) as mocked_uuid4:
            actual = record.apply_variant(mocked_variant, relative_start)
            assert mocked_uuid4.call_count == 0
        expected = []

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantEndsInNextLoci_skips(self, mocked_variant):
        name = "foo"
        comment = ""
        seq = "yyACG"
        record = Record(name, comment, seq)
        relative_start = 5
        mocked_variant.POS = 7
        mocked_variant.REF = "ACGTACGT"
        mocked_variant.ALT = ["abcdefgh"]

        with patch("uuid.uuid4", return_value=name) as mocked_uuid4:
            actual = record.apply_variant(mocked_variant, relative_start)
            assert mocked_uuid4.call_count == 0
        expected = []

        assert actual == expected
