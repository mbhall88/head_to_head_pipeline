from pathlib import Path


rule unicycler:
    input:
         illumina1  = rules.trim_illumina.output.forward_paired,
         illumina2  = rules.trim_illumina.output.reverse_paired,
         long_reads = mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
          assembly       = outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.unicycler.{technology}.fasta",
          assembly_graph = outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.unicycler.{technology}.gfa",
    threads: 16
    resources:
             mem_mb = lambda wildcards, attempt: 16000 * attempt
    singularity: containers["unicycler"]
    params:
        outdir = lambda wildcards, output: Path(output.assembly).parent,
        verbosity = 2,
        keep = 0,  # keep only assembly, graph, and log
    shell:
        """
        unicycler --short1 {input.illumina1} \
            --short2 {input.illumina2} \
            --long {input.long_reads} \
            --out {params.outdir} \
            --verbosity {params.verbosity} \
            --threads {threads} 
        mv {params.outdir}/assembly.fasta {output.assembly}
        mv {params.outdir}/assembly.gfa {output.assembly_graph}
        """
