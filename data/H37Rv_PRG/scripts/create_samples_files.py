with open(snakemake.output.samples_file, "w") as outstream, open(
    snakemake.input.assignments
) as instream:
    for row in map(str.rstrip, instream):
        fields = row.split(",")
        sample = fields[0]
        lin = fields[1]
        if (
            snakemake.wildcards.lineage == "rare"
            and lin.lower() in snakemake.params.rare_lineages
        ) or lin == snakemake.wildcards.lineage:
            print(sample, file=outstream)
