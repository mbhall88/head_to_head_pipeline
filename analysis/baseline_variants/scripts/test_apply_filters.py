from unittest.mock import patch

from apply_filters import *
from pytest import raises


class TestGetDepth:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noDepthTag_returnsDefault(self, mocked_variant):
        default_depth = 0
        mocked_variant.INFO = dict()

        actual = get_depth(mocked_variant, default=default_depth)
        expected = default_depth

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_depthTagPresent_returnsDepth(self, mocked_variant):
        dp = 6
        mocked_variant.INFO = {str(Tags.Depth): dp}

        actual = get_depth(mocked_variant)
        expected = dp

        assert actual == expected


class TestGetStrandDepths:
    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_noStrandDepthsTag_returnsDefault(self, mocked_variant):
        default_strand_depths = StrandDepths()
        mocked_variant.INFO = dict()

        actual = get_strand_depths(mocked_variant, default=default_strand_depths)
        expected = default_strand_depths

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandDepthTagPresent_returnsDepth(self, mocked_variant):
        strand_depth = StrandDepths(4, 5, 2, 8)
        mocked_variant.INFO = {str(Tags.StrandDepth): strand_depth.to_tuple()}

        actual = get_strand_depths(mocked_variant)
        expected = strand_depth

        assert actual == expected


class TestFilterStatusStr:
    def test_noFilters_returnsPass(self):
        status = FilterStatus()

        actual = str(status)
        expected = Tags.Pass.value

        assert actual == expected

    def test_allFilters_returnsAllInStr(self):
        delim = ";"
        status = FilterStatus(
            low_depth=True,
            high_depth=True,
            low_qual=True,
            strand_bias=True,
            delim=delim,
        )

        actual = str(status)
        expected = delim.join(
            map(str, [Tags.LowDepth, Tags.HighDepth, Tags.LowQual, Tags.StrandBias])
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
    def test_lowDepthOnAndVarHasLowDepth(self, mocked_variant):
        expected_depth = 100
        min_depth = 0.5
        assessor = Filter(expected_depth=expected_depth, min_depth=min_depth)
        variant_depth = 5

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_depth=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowDepthOnAndVarHasFineDepth(self, mocked_variant):
        expected_depth = 100
        min_depth = 0.5
        assessor = Filter(expected_depth=expected_depth, min_depth=min_depth)
        variant_depth = 90

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowDepthOnAndVarHasMinDepth(self, mocked_variant):
        expected_depth = 100
        min_depth = 0.5
        assessor = Filter(expected_depth=expected_depth, min_depth=min_depth)
        variant_depth = 50

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxDepthOnAndVarHasHighDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        assessor = Filter(expected_depth=expected_depth, max_depth=max_depth)
        variant_depth = 201

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_depth=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxDepthOnAndVarHasOkDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        assessor = Filter(expected_depth=expected_depth, max_depth=max_depth)
        variant_depth = expected_depth

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_maxDepthOnAndVarHasMaxDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        assessor = Filter(expected_depth=expected_depth, max_depth=max_depth)
        variant_depth = 200

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothDepthFiltersOnAndVarHasExpectedDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        min_depth = 0.5
        assessor = Filter(
            expected_depth=expected_depth, max_depth=max_depth, min_depth=min_depth
        )
        variant_depth = expected_depth

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(Tags.Pass)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothDepthFiltersOnAndMaxLowerThanMin_raisesError(self, mocked_variant):
        expected_depth = 100
        max_depth = 0.2
        min_depth = 0.5

        with raises(ValueError) as err:
            Filter(
                expected_depth=expected_depth, max_depth=max_depth, min_depth=min_depth
            )
            assert "Minimum depth is more than maximum depth" in err

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothDepthFiltersOnAndVarHasLowDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        min_depth = 0.5
        assessor = Filter(
            expected_depth=expected_depth, max_depth=max_depth, min_depth=min_depth
        )
        variant_depth = 2

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_depth=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_bothDepthFiltersOnAndVarHasHighDepth(self, mocked_variant):
        expected_depth = 100
        max_depth = 2.0
        min_depth = 0.5
        assessor = Filter(
            expected_depth=expected_depth, max_depth=max_depth, min_depth=min_depth
        )
        variant_depth = 250

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_depth=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualOnAndVarHasHighQual(self, mocked_variant):
        mocked_variant.QUAL = 50
        min_qual = 20
        assessor = Filter(min_qual=min_qual)

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualOnAndVarHasLowQual(self, mocked_variant):
        mocked_variant.QUAL = 50
        min_qual = 200
        assessor = Filter(min_qual=min_qual)

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_qual=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualOffAndVarHasLowQual(self, mocked_variant):
        mocked_variant.QUAL = 1
        assessor = Filter()

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualOnAndVarHasMinQual(self, mocked_variant):
        min_qual = 200
        mocked_variant.QUAL = min_qual
        assessor = Filter(min_qual=min_qual)

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualOnAndVarHasMinQualMinusOne(self, mocked_variant):
        min_qual = 200
        mocked_variant.QUAL = min_qual - 1
        assessor = Filter(min_qual=min_qual)

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_qual=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowQualLowDepthOnAndVarFailsBoth(self, mocked_variant):
        min_qual = 200
        mocked_variant.QUAL = 3
        expected_depth = 100
        min_depth = 0.5
        assessor = Filter(
            expected_depth=expected_depth, min_depth=min_depth, min_qual=min_qual
        )
        variant_depth = 4

        with patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_qual=True, low_depth=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasNoBias(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(7, 7, 0, 0)
        mocked_variant.genotypes = [[0]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasBias(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(1, 7, 0, 0)
        mocked_variant.genotypes = [[0]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasNoRefDepth(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 0, 0, 0)
        mocked_variant.genotypes = [[0]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasNoAltDepth(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 0, 0, 0)
        mocked_variant.genotypes = [[1]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarHasBiasOnAlt(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(1, 1, 1, 7)
        mocked_variant.genotypes = [[1]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarIsNull(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 0, 0, 0)
        mocked_variant.genotypes = [[-1]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarBiasIsOnLimit(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 0, 5, 15)
        mocked_variant.genotypes = [[1]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarBiasIsOneBelowLimit(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 0, 24, 76)
        mocked_variant.genotypes = [[1]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_strandBiasOnAndVarIsHet(self, mocked_variant):
        bias = 25
        assessor = Filter(min_strand_bias=bias)
        strand_depths = StrandDepths(0, 5, 29, 76)
        mocked_variant.genotypes = [[1, 0]]

        with patch("apply_filters.get_strand_depths", return_value=strand_depths):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(strand_bias=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_allFiltersOnAndVarPassesAll(self, mocked_variant):
        bias = 25
        strand_depths = StrandDepths(5, 5, 0, 0)
        mocked_variant.genotypes = [[0]]
        min_qual = 20
        mocked_variant.QUAL = 30
        expected_depth = 10
        min_depth = 0.1
        max_depth = 2.0
        assessor = Filter(
            expected_depth=expected_depth,
            min_depth=min_depth,
            max_depth=max_depth,
            min_qual=min_qual,
            min_strand_bias=bias,
        )
        variant_depth = 10

        with patch(
            "apply_filters.get_strand_depths", return_value=strand_depths
        ), patch("apply_filters.get_depth", return_value=variant_depth):
            actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus())

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowBaseQualBiasAndVarHasLowBaseQualBias(self, mocked_variant):
        min_bqb = 0.5
        assessor = Filter(min_bqb=min_bqb)
        bqb = 0.1
        mocked_variant.INFO = {str(Tags.BaseQualBias): bqb}

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_bqb=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowMappingQualBiasAndVarHasLowMappingQualBias(self, mocked_variant):
        min_mqb = 0.5
        assessor = Filter(min_mqb=min_mqb)
        mqb = 0.1
        mocked_variant.INFO = {str(Tags.MapQualBias): mqb}

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_mqb=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowReadPosBiasAndVarHasLowReadPosBias(self, mocked_variant):
        min_rpb = 0.5
        assessor = Filter(min_rpb=min_rpb)
        rpb = 0.1
        mocked_variant.INFO = {str(Tags.ReadPosBias): rpb}

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_rpb=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_highSegBiasAndVarHasHighSegBias(self, mocked_variant):
        max_sgb = -0.5
        assessor = Filter(max_sgb=max_sgb)
        sgb = -0.4
        mocked_variant.INFO = {str(Tags.SegregationBias): sgb}

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(high_sgb=True))

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_lowVarDistBiasAndVarHasLowVarDistBias(self, mocked_variant):
        min_vdb = 0.5
        assessor = Filter(min_vdb=min_vdb)
        vdb = 0.1
        mocked_variant.INFO = {str(Tags.VariantDistanceBias): vdb}

        actual = assessor.filter_status(mocked_variant)
        expected = str(FilterStatus(low_vdb=True))

        assert actual == expected
