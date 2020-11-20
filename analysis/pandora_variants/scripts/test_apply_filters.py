from unittest.mock import patch, call

from apply_filters import *
from pytest import raises


class TestGetCovg:
    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_covgTagPresent_returnsCovg(self, mocked_variant, mocked_strand):
        fwd_covg = 6
        rev_covg = 1
        mocked_strand.return_value = Strand(fwd_covg, rev_covg)

        actual = get_covg(mocked_variant)
        expected = fwd_covg + rev_covg

        assert actual == expected


@patch("cyvcf2.Variant", autospec=True, create=True)
def test_get_gt_conf(mocked_variant):
    expected = 20.7
    mocked_variant.format.return_value = [[expected]]

    actual = get_gt_conf(mocked_variant)

    assert actual == expected

    calls = [
        call(Tags.GtypeConf.value),
    ]
    mocked_variant.format.assert_has_calls(calls, any_order=False)


@patch("cyvcf2.Variant", autospec=True, create=True)
def test_get_gaps(mocked_variant):
    mocked_variant.format.return_value = [[0, 0.5]]
    mocked_variant.genotypes = [[1]]

    actual = get_gaps(mocked_variant)
    expected = 0.5

    assert actual == expected

    calls = [
        call(Tags.Gaps.value),
    ]
    mocked_variant.format.assert_has_calls(calls, any_order=False)


class TestGetStrandCovgs:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noCovgTag_raiseError(self, mocked_variant):
        mocked_variant.format.side_effect = KeyError
        mocked_variant.genotypes = [[0]]

        with raises(KeyError):
            strand = Strand.from_variant(mocked_variant)

        calls = [
            call(Tags.FwdCovg.value),
        ]
        mocked_variant.format.assert_has_calls(calls, any_order=False)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_covgTagsPresent_returnsCovgs(self, mocked_variant):
        fwd_covg = 6
        rev_covg = 2
        mocked_variant.format.side_effect = [[[fwd_covg]], [[rev_covg]]]
        mocked_variant.genotypes = [[0]]
        strand = Strand.from_variant(mocked_variant)

        actual = strand.covgs
        expected = (fwd_covg, rev_covg)

        assert actual == expected

        calls = [
            call(Tags.FwdCovg.value),
            call(Tags.RevCovg.value),
        ]
        mocked_variant.format.assert_has_calls(calls, any_order=False)


class TestFilterStatusStr:
    def test_noFilters_returnsPass(self):
        status = FilterStatus()

        actual = str(status)
        expected = Tags.Pass.value

        assert actual == expected

    def test_allFilters_returnsAllInStr(self):
        delim = ";"
        status = FilterStatus(
            low_covg=True,
            high_covg=True,
            high_gaps=True,
            strand_bias=True,
            low_gt_conf=True,
            delim=delim,
        )

        actual = str(status)
        expected = delim.join(
            map(
                str,
                [
                    Tags.LowCovg,
                    Tags.HighCovg,
                    Tags.LowGtConf,
                    Tags.StrandBias,
                    Tags.HighGaps,
                ],
            )
        )

        assert actual == expected

    def test_oneFilter_returnsOneWithoutDelim(self):
        delim = ";"
        status = FilterStatus(strand_bias=True, delim=delim,)

        actual = str(status)
        expected = str(Tags.StrandBias)

        assert actual == expected


class TestFilterFilterStatus:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noFilters_returnsPass(self, mocked_variant):
        assessor = Filter()

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowCovgOnAndVarHasLowCovg(self, mocked_variant):
        min_covg = 10
        assessor = Filter(min_covg=min_covg)
        variant_covg = 9

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_covg=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowCovgOnAndVarHasGoodCovg(self, mocked_variant):
        min_covg = 9
        assessor = Filter(min_covg=min_covg)
        variant_covg = 12

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowCovgOnAndVarHasMinCovg(self, mocked_variant):
        min_covg = 9
        assessor = Filter(min_covg=min_covg)
        variant_covg = 9

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxCovgOnAndVarHasHighCovg(self, mocked_variant):
        max_covg = 20
        assessor = Filter(max_covg=max_covg)
        variant_covg = 201

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_covg=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxCovgOnAndVarHasOkCovg(self, mocked_variant):
        max_covg = 20
        assessor = Filter(max_covg=max_covg)
        variant_covg = 2

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxCovgOnAndVarHasMaxCovg(self, mocked_variant):
        max_covg = 20
        assessor = Filter(max_covg=max_covg)
        variant_covg = 20

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothCovgFiltersOnAndVarHasExpectedCovg(self, mocked_variant):
        min_covg = 5
        max_covg = 20
        assessor = Filter(max_covg=max_covg, min_covg=min_covg)
        variant_covg = 15

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothCovgFiltersOnAndMaxLowerThanMin_raisesError(self, mocked_variant):
        expected_covg = 100
        max_covg = 0.2
        min_covg = 0.5

        with raises(ValueError) as err:
            Filter(max_covg=max_covg, min_covg=min_covg)
            assert "Minimum covg is more than maximum covg" in err

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothCovgFiltersOnAndVarHasLowCovg(self, mocked_variant):
        min_covg = 5
        max_covg = 20
        assessor = Filter(max_covg=max_covg, min_covg=min_covg)
        variant_covg = 1

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_covg=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothCovgFiltersOnAndVarHasHighCovg(self, mocked_variant):
        min_covg = 5
        max_covg = 20
        assessor = Filter(max_covg=max_covg, min_covg=min_covg)
        variant_covg = 100

        with patch("apply_filters.get_covg", return_value=variant_covg):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_covg=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowGtConfOnAndVarHasHighGtConf(self, mocked_variant):
        min_gt_conf = 1.1
        assessor = Filter(min_gt_conf=min_gt_conf)
        variant_gt_conf = 2.2

        with patch("apply_filters.get_gt_conf", return_value=variant_gt_conf):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowGtConfOnAndVarHasLowGtConf(self, mocked_variant):
        min_gt_conf = 1.1
        assessor = Filter(min_gt_conf=min_gt_conf)
        variant_gt_conf = 0.5

        with patch("apply_filters.get_gt_conf", return_value=variant_gt_conf):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_gt_conf=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowGtConfOnAndVarHasMinGtConf(self, mocked_variant):
        min_gt_conf = 1.1
        assessor = Filter(min_gt_conf=min_gt_conf)
        variant_gt_conf = 1.1

        with patch("apply_filters.get_gt_conf", return_value=variant_gt_conf):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowGtConfOnAndVarHasMinGtConfMinusOne(self, mocked_variant):
        min_gt_conf = 1.1
        assessor = Filter(min_gt_conf=min_gt_conf)
        variant_gt_conf = 1.09

        with patch("apply_filters.get_gt_conf", return_value=variant_gt_conf):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_gt_conf=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowGtConfLowCovgOnAndVarFailsBoth(self, mocked_variant):
        min_gt_conf = 1.1
        min_covg = 9
        assessor = Filter(min_covg=min_covg, min_gt_conf=min_gt_conf)
        variant_covg = 4
        variant_gt_conf = 0.8

        with patch("apply_filters.get_covg", return_value=variant_covg), patch(
            "apply_filters.get_gt_conf", return_value=variant_gt_conf
        ):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_gt_conf=True, low_covg=True))

        assert actual == expected

    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasNoBias(self, mocked_variant, mocked_strand):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_covgs = Strand(5, 5)
        mocked_strand.return_value = strand_covgs

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasBias(self, mocked_variant, mocked_strand):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_covgs = Strand(5, 1)
        mocked_strand.return_value = strand_covgs

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasNoCovg(self, mocked_variant, mocked_strand):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_covgs = Strand(0, 0)
        mocked_strand.return_value = strand_covgs

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarIsNull_usesRef(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        mocked_variant.genotypes = [[-1]]
        fwd_covg = 6
        rev_covg = 1
        mocked_variant.format.side_effect = [[[fwd_covg, 0]], [[rev_covg, 0]]]

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarBiasIsOnLimit(self, mocked_variant, mocked_strand):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_covgs = Strand(15, 5)
        mocked_strand.return_value = strand_covgs

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch.object(Strand, Strand.from_variant.__name__)
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarBiasIsOneBelowLimit(self, mocked_variant, mocked_strand):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_covgs = Strand(24, 76)
        mocked_strand.return_value = strand_covgs

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarIsHet(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        mocked_variant.genotypes = [[0, 1]]

        with raises(NotImplementedError):
            assessor.filter_status(mocked_variant)

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxGapsOnAndVarHasLowGaps(self, mocked_variant):
        max_gaps = 0.5
        assessor = Filter(max_gaps=max_gaps)
        variant_gaps = 0.3

        with patch("apply_filters.get_gaps", return_value=variant_gaps):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxGapsOnAndVarHasHighGaps(self, mocked_variant):
        max_gaps = 0.5
        assessor = Filter(max_gaps=max_gaps)
        variant_gaps = 0.6

        with patch("apply_filters.get_gaps", return_value=variant_gaps):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_gaps=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxIndelOnAndIndelIsLong(self, mocked_variant):
        max_indel = 2
        assessor = Filter(max_indel=max_indel)
        mocked_variant.REF = "AAAA"
        mocked_variant.ALT = ["ATAA", "A"]
        mocked_variant.genotypes = [[2]]

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(long_indel=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxIndelOnAndIndelIsSameAsMax(self, mocked_variant):
        max_indel = 3
        assessor = Filter(max_indel=max_indel)
        mocked_variant.REF = "AAAA"
        mocked_variant.ALT = ["ATAA", "A"]
        mocked_variant.genotypes = [[2]]

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(long_indel=False))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxIndelOnAndIndelIsBelowMax(self, mocked_variant):
        max_indel = 1
        assessor = Filter(max_indel=max_indel)
        mocked_variant.REF = "AAAA"
        mocked_variant.ALT = ["ATAA", "A"]
        mocked_variant.genotypes = [[1]]

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(long_indel=False))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxIndelOffAndIndelIsLong(self, mocked_variant):
        max_indel = None
        assessor = Filter(max_indel=max_indel)
        mocked_variant.REF = "AA"
        mocked_variant.ALT = ["ATACGA", "A"]
        mocked_variant.genotypes = [[1]]

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(long_indel=False))

        assert actual == expected
