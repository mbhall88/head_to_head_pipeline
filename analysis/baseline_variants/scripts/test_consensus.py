from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner
from consensus import *


class TestBed:
    def test_invalidBedFormat_raisesError(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "foo bar"
        bedfile.write_text(content)

        with pytest.raises(ValueError):
            Bed(bedfile)

    def test_variantNotInBedInterval(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 8

        assert mocked_variant not in bed

    def test_variantPosInBedButDiffChrom(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr2"
        mocked_variant.POS = 3

        assert mocked_variant not in bed

    def test_variantPosAndChromInBed(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 4

        assert mocked_variant in bed

    def test_passSomethingNotVariant_raisesError(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)

        with pytest.raises(AttributeError):
            "foo" in bed

    def test_bedStartIsInclusive(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 3

        assert mocked_variant in bed

    def test_bedStartIsInclusiveAlsoWithZeroBased(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t3"
        bedfile.write_text(content)
        bed = Bed(bedfile, zero_based=False)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 2

        assert mocked_variant in bed

    def test_bedEndIsNotInclusive(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t5"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 6

        assert mocked_variant not in bed

    def test_bedStartIsNotInclusiveAlsoWithZeroBased(self, tmp_path: Path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t3"
        bedfile.write_text(content)
        bed = Bed(bedfile, zero_based=False)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 3

        assert mocked_variant not in bed


class TestClassifier:
    def test_variantInMaskAndIgnoreMask(self, tmp_path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t3"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 3
        mocked_variant.REF = "A"
        classifier = Classifier(mask=bed, ignore_mask=True)

        actual = classifier.classify(mocked_variant)
        expected = N

        assert actual == expected

    def test_variantInMaskAndNotIgnoreMask(self, tmp_path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t3"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        ref = "A"
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 3
        mocked_variant.genotypes = [[0]]
        mocked_variant.REF = ref
        classifier = Classifier(mask=bed, ignore_mask=False)

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    def test_variantNotInMaskAndIgnoreMask(self, tmp_path):
        bedfile = tmp_path / "tmp.bed"
        content = "chr1\t2\t3"
        bedfile.write_text(content)
        bed = Bed(bedfile)
        ref = "A"
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.CHROM = "chr1"
        mocked_variant.POS = 9
        mocked_variant.genotypes = [[0]]
        mocked_variant.REF = ref
        classifier = Classifier(mask=bed, ignore_mask=True)

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantFailsFilterAndIgnoreFilter(self, mocked_variant):
        classifier = Classifier(ignore_filter=True)
        mocked_variant.FILTER = "FAIL"
        mocked_variant.REF = "A"

        actual = classifier.classify(mocked_variant)
        expected = N

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantFailsFilterAndNotIgnoreFilter(self, mocked_variant):
        classifier = Classifier(ignore_filter=False)
        ref = "A"
        mocked_variant.FILTER = "FAIL"
        mocked_variant.genotypes = [[0]]
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantPassFilterAndIgnoreFilter(self, mocked_variant):
        classifier = Classifier(ignore_filter=False)
        ref = "A"
        mocked_variant.FILTER = None
        mocked_variant.genotypes = [[0]]
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsNullAndIgnoreNull(self, mocked_variant):
        classifier = Classifier(ignore_null=True)
        mocked_variant.genotypes = [[-1]]
        mocked_variant.REF = "A"

        actual = classifier.classify(mocked_variant)
        expected = N

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsNullAndNotIgnoreNull(self, mocked_variant):
        classifier = Classifier(ignore_null=False)
        ref = "T"
        mocked_variant.genotypes = [[-1]]
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsHetAndDefaultIsNone(self, mocked_variant):
        classifier = Classifier(het_default="none")
        mocked_variant.genotypes = [[1, 0]]
        mocked_variant.REF = "A"

        actual = classifier.classify(mocked_variant)
        expected = N

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsHetAndDefaultIsRef(self, mocked_variant):
        classifier = Classifier(het_default="ref")
        ref = "C"
        mocked_variant.genotypes = [[1, 0]]
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    def test_variantIsHetAndDefaultIsUnknown_rasiesError(self):
        with pytest.raises(UnknownDefaultHet):
            classifier = Classifier(het_default="unknown")

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsRefAltHetAndDefaultIsAlt(self, mocked_variant):
        classifier = Classifier(het_default="alt")
        ref = "C"
        alt = ["G"]
        mocked_variant.genotypes = [[1, 0]]
        mocked_variant.REF = ref
        mocked_variant.ALT = alt

        actual = classifier.classify(mocked_variant)
        expected = alt[0]

        assert actual == expected

    def test_variantIsAltAltHetAndDefaultIsAlt(self, caplog):
        classifier = Classifier(het_default="alt")
        alt = ["G", "T"]
        mocked_variant = patch("cyvcf2.Variant", autospec=True, create=True)
        mocked_variant.genotypes = [[1, 2]]
        mocked_variant.ALT = alt
        mocked_variant.REF = "A"

        actual = classifier.classify(mocked_variant)
        expected = alt[1]

        assert actual == expected

        assert "Using ALT with highest index" in caplog.text

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsAltAltHetAndDefaultIsRef(self, mocked_variant):
        classifier = Classifier(het_default="ref")
        ref = "A"
        alt = ["G", "T"]
        mocked_variant.genotypes = [[1, 2]]
        mocked_variant.ALT = alt
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsRefEverythingElseDefault(self, mocked_variant):
        classifier = Classifier()
        ref = "A"
        mocked_variant.genotypes = [[0, 0]]
        mocked_variant.REF = ref

        actual = classifier.classify(mocked_variant)
        expected = ref

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsMNPButFailsFilter(self, mocked_variant):
        classifier = Classifier(ignore_filter=True)
        ref = "AT"
        mocked_variant.genotypes = [[0, 0]]
        mocked_variant.REF = ref
        mocked_variant.FILTER = "FAIL"

        actual = classifier.classify(mocked_variant)
        expected = N * len(ref)

        assert actual == expected

    @patch("cyvcf2.Variant", autospec=True, create=True)
    def test_variantIsAltEverythingElseDefault(self, mocked_variant):
        classifier = Classifier()
        ref = "A"
        alt = ["G"]
        mocked_variant.genotypes = [[1, 1]]
        mocked_variant.REF = ref
        mocked_variant.ALT = alt

        actual = classifier.classify(mocked_variant)
        expected = alt[0]

        assert actual == expected


def create_variant(
    pos, ref, alt, gt, chrom="s1", _id=".", qual=285, filt="PASS", info=".", fmt="GT"
) -> str:
    return "\t".join(
        [chrom, str(pos), _id, ref, ",".join(alt), str(qual), filt, info, fmt, gt]
    )


class TestMain:
    @pytest.fixture
    def in_vcf(self, tmp_path):
        header = "\t".join(
            [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                "sample1",
            ]
        )

        vcf_content = f"""##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FAIL,Description="All filters failed">
##contig=<ID=s1,length=23>
#{header}
{create_variant(pos=3, ref="T", alt=["x"], gt="1", filt="FAIL")}
{create_variant(pos=5, ref="G", alt=["."], gt=".")}
{create_variant(pos=9, ref="C", alt=["A"], gt="0/1")}
{create_variant(pos=10, ref="G", alt=["C", "A"], gt="1/2")}
{create_variant(pos=12, ref="T", alt=["x"], gt="1")}
{create_variant(pos=14, ref="G", alt=["C"], gt="0")}
{create_variant(pos=18, ref="C", alt=["x"], gt="1")}
{create_variant(pos=23, ref="C", alt=["x"], gt="1/.")}
"""
        vcfpath = tmp_path / "tmp.vcf"
        vcfpath.write_text(vcf_content)
        return str(vcfpath)

    @pytest.fixture
    def mask(self, tmp_path):
        content = "\n".join(["s1\t17\t18", "s1\t20\t22"])
        bedpath = tmp_path / "tmp.bed"
        bedpath.write_text(content)
        return str(bedpath)

    @pytest.fixture
    def ref_fasta(self, tmp_path):
        content = ">s1\nAATCGAACCGTTTGGAGCAGTTC"
        fastapath = tmp_path / "tmp.fa"
        fastapath.write_text(content)
        return str(fastapath)

    def test_ignoreAllAndHets(self, in_vcf, mask, ref_fasta):
        ignore = "all"
        het_default = "none"
        runner = CliRunner()
        args = [
            "-i",
            in_vcf,
            "-f",
            ref_fasta,
            "-m",
            mask,
            "-H",
            het_default,
            "-I",
            ignore,
        ]
        result = runner.invoke(main, args,)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nNNNNNNNNNNNxNGNNNNNNNNx\n"

        assert actual == expected

    def test_ignoreAllAndRefHetDefault(self, in_vcf, mask, ref_fasta):
        ignore = "all"
        het_default = "ref"
        runner = CliRunner()
        args = [
            "-i",
            in_vcf,
            "-f",
            ref_fasta,
            "-m",
            mask,
            "-H",
            het_default,
            "-I",
            ignore,
        ]
        result = runner.invoke(main, args,)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nNNNNNNNNCGNxNGNNNNNNNNx\n"

        assert actual == expected

    def test_ignoreAllAndAltHetDefault(self, in_vcf, mask, ref_fasta):
        ignore = "all"
        het_default = "alt"
        runner = CliRunner()
        args = [
            "-i",
            in_vcf,
            "-f",
            ref_fasta,
            "-m",
            mask,
            "-H",
            het_default,
            "-I",
            ignore,
        ]
        result = runner.invoke(main, args,)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nNNNNNNNNAANxNGNNNNNNNNx\n"

        assert actual == expected

    def test_ignoreNullAndFilterAndIgnoreHet(self, in_vcf, mask, ref_fasta):
        ignore = ["filter", "null"]
        het_default = "none"
        runner = CliRunner()
        args = ["-i", in_vcf, "-f", ref_fasta, "-m", mask, "-H", het_default]
        for opt in ignore:
            args.extend(["-I", opt])
        result = runner.invoke(main, args)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nAANCNAACNNTxTGGAGxAGTTx\n"

        assert actual == expected

    def test_ignoreNoneAndRefHetDefault(self, in_vcf, mask, ref_fasta):
        ignore = ["none"]
        het_default = "ref"
        runner = CliRunner()
        args = ["-i", in_vcf, "-f", ref_fasta, "-m", mask, "-H", het_default]
        for opt in ignore:
            args.extend(["-I", opt])
        result = runner.invoke(main, args)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nAAxCGAACCGTxTGGAGxAGTTx\n"

        assert actual == expected

    def test_ignoreMissingAndMaskAltHetDefault(self, in_vcf, mask, ref_fasta):
        ignore = ["missing", "mask"]
        het_default = "alt"
        runner = CliRunner()
        args = ["-i", in_vcf, "-f", ref_fasta, "-m", mask, "-H", het_default]
        for opt in ignore:
            args.extend(["-I", opt])
        result = runner.invoke(main, args)

        assert result.exit_code == 0

        actual = result.output
        expected = ">sample1\nNNxNGNNNAANxNGNNNNNNNNx\n"

        assert actual == expected

    def test_ignoreNoneAndMissing_raisesError(self, in_vcf, mask, ref_fasta):
        ignore = ["missing", "none"]
        het_default = "alt"
        runner = CliRunner()
        args = ["-i", in_vcf, "-f", ref_fasta, "-m", mask, "-H", het_default]
        for opt in ignore:
            args.extend(["-I", opt])
        result = runner.invoke(main, args)

        assert result.exit_code != 0
