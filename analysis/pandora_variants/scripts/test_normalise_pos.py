from pathlib import Path
from click.testing import CliRunner

from normalise_pos import *


def create_vcf(chrom: str, pos: int, ref: str) -> Path:
    tmpfile = Path("tmp.vcf")
    header = "##fileformat=VCFv4.3\n"
    header += "\t".join(
        "#CHROM  POS  ID  REF  ALT QUAL  FILTER  INFO  FORMAT  sample".split()
    )
    var = "\t".join(f"{chrom}  {pos}  .   {ref} . .   .   . GT  0".split())
    tmpfile.write_text(f"{header}\n{var}\n")
    return tmpfile


def create_ref(header: str, seq: str) -> Path:
    tmpfile = Path("tmp.fa")
    tmpfile.write_text(f">{header}\n{seq}\n")
    return tmpfile


def loci_info() -> Path:
    header = "filename,type,start,end,name,contig"
    row = "filename,type,10,13,loci1,contig"
    tmpfile = Path("tmp.loci_info.csv")
    tmpfile.write_text(f"{header}\n{row}\n")
    return tmpfile


def test_loci_in_vcf_does_not_exist_in_info():
    locifile = loci_info()
    chrom = "nonexistent"
    pos = 1
    ref = "A"
    vcffile = create_vcf(chrom=chrom, pos=pos, ref=ref)
    params = ["-i", str(vcffile), "-l", str(locifile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()

    assert result.exit_code != 0
    expected_exception = KeyError(f"{chrom} is not in the loci information file")
    assert type(result.exception) == type(expected_exception)
    assert str(result.exception) == str(expected_exception)


def test_pos_outside_expected_length():
    locifile = loci_info()
    chrom = "loci1"
    pos = 4
    ref = "A"
    vcffile = create_vcf(chrom=chrom, pos=pos, ref=ref)
    params = ["-i", str(vcffile), "-l", str(locifile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()

    assert result.exit_code != 0
    expected_exception = IndexError(
        f"Position {pos} for loci {chrom} is outside the expected length of the loci"
    )
    assert type(result.exception) == type(expected_exception)
    assert str(result.exception) == str(expected_exception)


def test_original_genome_seq_doesnt_match_at_only_pos():
    locifile = loci_info()
    chrom = "loci1"
    pos = 3
    ref = "T"
    vcffile = create_vcf(chrom=chrom, pos=pos, ref=ref)
    seq = "AAA"
    ref_header = "contig"
    reffile = create_ref(header=ref_header, seq=seq)
    params = ["-i", str(vcffile), "-l", str(locifile), "-r", str(reffile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()
    reffile.unlink()
    print(result.exception)
    assert result.exit_code != 0
    expected_exception = ReferenceError(
        f"VCF REF {ref} does not match the expected reference "
        f"sequence A at loci {chrom} position {pos}"
    )
    assert type(result.exception) == type(expected_exception)
    assert str(result.exception) == str(expected_exception)


def test_original_genome_seq_doesnt_match_at_second_pos():
    locifile = loci_info()
    chrom = "loci1"
    pos = 2
    ref = "AT"
    vcffile = create_vcf(chrom=chrom, pos=pos, ref=ref)
    seq = "AAA"
    ref_header = "contig"
    reffile = create_ref(header=ref_header, seq="AAA")
    params = ["-i", str(vcffile), "-l", str(locifile), "-r", str(reffile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()
    reffile.unlink()

    assert result.exit_code != 0
    expected_exception = ReferenceError(
        f"VCF REF {ref} does not match the expected reference "
        f"sequence AA at loci {chrom} position {pos}"
    )
    assert type(result.exception) == type(expected_exception)
    assert str(result.exception) == str(expected_exception)


def test_loci_contig_not_in_original_genome():
    locifile = loci_info()
    chrom = "loci1"
    pos = 2
    ref = "AT"
    vcffile = create_vcf(chrom=chrom, pos=pos, ref=ref)
    seq = "AAA"
    ref_header = "foo"
    reffile = create_ref(header=ref_header, seq="AAA")
    params = ["-i", str(vcffile), "-l", str(locifile), "-r", str(reffile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()
    reffile.unlink()

    assert result.exit_code != 0
    expected_exception = IndexError(f"Contig contig not in the reference genome")
    assert type(result.exception) == type(expected_exception)
    assert str(result.exception) == str(expected_exception)


def test_all_fields_error_free():
    locifile = loci_info()
    chrom = "loci1"
    loci_pos = 2
    ref = "AA"
    vcffile = create_vcf(chrom=chrom, pos=loci_pos, ref=ref)
    seq = "AAA"
    ref_header = "contig"
    reffile = create_ref(header=ref_header, seq="AAA")
    params = ["-i", str(vcffile), "-l", str(locifile), "-r", str(reffile)]
    runner = CliRunner()
    result = runner.invoke(main, params)

    locifile.unlink()
    vcffile.unlink()
    reffile.unlink()

    assert result.exit_code == 0

    actual = result.output
    expected_chrom = "contig"
    expected_pos = 12
    expected_loci = chrom
    expected_loci_pos = loci_pos
    expected_ref = ref
    expected = f"""##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=LOCI,Number=1,Type=String,Description="CHROM (loci) in the original VCF">
##INFO=<ID=LPOS,Number=1,Type=Integer,Description="LOCI (see other INFO header) POS in the original VCF">
##contig=<ID={expected_chrom}>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
{expected_chrom}\t{expected_pos}\t.\t{expected_ref}\t.\t.\t.\tLOCI={expected_loci};LPOS={expected_loci_pos}\tGT\t0
"""

    assert actual == expected


def test_reset_contigs_in_header():
    header_col = "\t".join(
        "#CHROM  POS  ID  REF  ALT QUAL  FILTER  INFO  FORMAT  sample".split()
    )
    header = (
        "##fileformat=VCFv4.3\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "##contig=<ID=old>\n"
        f"{header_col}\n"
    )
    contigs = sorted({"foo", "bar"})

    actual = reset_contigs_in_header(header, contigs)
    expected = (
        "##fileformat=VCFv4.3\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "##contig=<ID=bar>\n"
        "##contig=<ID=foo>\n"
        f"{header_col}"
    )

    assert actual == expected
