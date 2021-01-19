import sys

sys.stderr = open(snakemake.log[0], "w")
from pathlib import Path
from typing import List
import uuid
from collections import defaultdict
from itertools import chain
import fileinput


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def extract_prg_name(filepath: Path) -> str:
    name = filepath.name
    return name.split(".")[0]


def concatenate(infiles: List[Path], outfile: Path):
    with outfile.open("w") as outstream, fileinput.input(files=infiles) as instream:
        for line in map(str.rstrip, instream):
            if line.startswith(">"):
                old_header = line[1:]
                new_header = f">{uuid.uuid4()} {old_header}"
                print(new_header, file=outstream)
            else:
                print(line, file=outstream)


input_dirs = [Path(p) for p in snakemake.input.denovo_dirs]
eprint(f"Aggregating de novo paths from {len(input_dirs)} directories")
outdir = Path(snakemake.output.outdir)
if not outdir.exists():
    outdir.mkdir()
input_files = chain.from_iterable(d.glob("*.fa") for d in input_dirs)

prg_name_to_paths = defaultdict(list)
for filepath in input_files:
    prg_name = extract_prg_name(filepath)
    prg_name_to_paths[prg_name].append(filepath)

eprint(f"Concatenating sequences for {len(prg_name_to_paths)} PRGs...")
for prg_name, paths in prg_name_to_paths.items():
    outfile = outdir / f"{prg_name}.fa"
    concatenate(paths, outfile)

eprint("Done!")
