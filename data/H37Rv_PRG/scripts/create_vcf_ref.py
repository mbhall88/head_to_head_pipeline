from pathlib import Path

outstream = open(snakemake.output.vcf_ref, "w")

with open(snakemake.input.loci_info) as csvfile:
    for row in map(str.rstrip, csvfile):
        fields = row.split(",")
        loci_fasta = Path(fields[0])
        outstream.write(loci_fasta.read_text())

outstream.close()
