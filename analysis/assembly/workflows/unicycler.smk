from pathlib import Path

rule unicycler:
    input:
        illumina1=rules.trim_illumina.output.forward_paired,
        illumina2=rules.trim_illumina.output.reverse_paired,
        long_reads=mada_dir / "{technology}" / "{sample}" / "{sample}.{technology}.fastq.gz",
    output:
        assembly=outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.unicycler.{technology}.fasta",
        assembly_graph=outdir / "{sample}" / "unicycler" / "{technology}" / "assembly.unicycler.{technology}.gfa",
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt
    singularity: containers["unicycler"]
    params:
        outdir=lambda wildcards, output: Path(output.assembly).parent,
        verbosity=2,
        keep=0,  # keep only assembly, graph, and log
        original_assembly_graph_name="assembly.gfa",
        original_assembly_name="assembly.fasta",
    shell:
        """
        unicycler --short1 {input.illumina1} \
            --short2 {input.illumina2} \
            --long {input.long_reads} \
            --out {params.outdir} \
            --verbosity {params.verbosity} \
            --threads {threads} 
        sleep 5
        cp {params.outdir}/{params.original_assembly_name} {output.assembly}
        cp {params.outdir}/{params.original_assembly_graph_name} {output.assembly_graph}
        """
