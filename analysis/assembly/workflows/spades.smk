from pathlib import Path


rule spades:
    input:
        illumina1 = rules.trim_illumina.output.forward_paired,
        illumina2 = rules.trim_illumina.output.reverse_paired,
        pacbio    = pacbio_dir / "{sample}" / "{sample}.pacbio.fastq.gz",
        nanopore  = ont_dir / "{sample}" / "{sample}.nanopore.fastq.gz"
    output:
        assembly       = outdir / "{sample}" / "spades" / "assembly.spades.{technology}.fasta",
        assembly_graph = outdir / "{sample}" / "spades" / "assembly_graph.spades.{technology}.gfa",
    shadow: "shallow"
    threads: 32
    resources:
        mem_mb = lambda wildcards, attempt: 64000 + (16000 * (attempt - 1)),
    params:
        mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1000),
        outdir = lambda wildcards, output: Path(output.assembly).parent
    singularity: containers["spades"]
    shell:
        """
        spades.py \
           -o {params.outdir} \
           --pe1-1 {input.illumina1} \
           --pe1-2 {input.illumina2} \
           --s1 {input.pacbio} \
           --nanopore {input.nanopore} \
           --threads {threads} \
           --memory {params.mem_gb} \
           --isolate
        cat {params.outdir}/spades.log
        mv {params.outdir}/scaffolds.fasta {output.assembly}
        mv {params.outdir}/assembly_graph.gfa {output.assembly_graph}
        """

