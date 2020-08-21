from pathlib import Path

outstream = open(snakemake.output.vcf_ref, "w")
loci_dir = Path(snakemake.input.loci_info).parent

with open(snakemake.input.loci_info) as csvfile:
    _ = next(csvfile)  # skip header
    for row in map(str.rstrip, csvfile):
        fields = row.split(",")
        loci_fasta = loci_dir / fields[0]
        print(loci_fasta.read_text(), file=outstream)

outstream.close()
