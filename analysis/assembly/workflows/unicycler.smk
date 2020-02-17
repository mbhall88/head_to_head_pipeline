rule unicycler:
    input:
         illumina1          = outdir / "{sample}" / "trimmed" / "{sample}.R1.trimmed.fastq.gz",
         illumina2          = outdir / "{sample}" / "trimmed" / "{sample}.R2.trimmed.fastq.gz",
         long_reads         = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
          assembly = outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.fasta",
          assembly_graph = outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.gfa",
    threads: 16
    resources:
             mem_mb = lambda wildcards, attempt: 16000 * attempt
    singularity: containers["unicycler"]
    params:
          verbosity = 2,
    shell:
        """
        outdir=$(dirname {output.assembly})
        unicycler --short1 {input.illumina1} \
            --short2 {input.illumina2} \
            --long {input.long_reads} \
            --out $outdir \
            --verbosity {params.verbosity} \
            --threads {threads} 
        """

# todo: add rule to annotate
# todo: add rules to polish
# todo: add rule to analyse pileup
