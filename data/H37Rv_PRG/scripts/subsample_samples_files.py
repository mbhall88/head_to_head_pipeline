import random
import sys
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

random.seed(snakemake.params.seed)
samples = list(
    filter(None, Path(snakemake.input.samples_file).read_text().splitlines())
)
threshold = snakemake.params.prg_names[snakemake.wildcards.prg_name]
k = min(int(threshold), len(samples))
selections = random.sample(samples, k=k)
content = "\n".join(selections)
Path(snakemake.output.samples_file).write_text(content)
