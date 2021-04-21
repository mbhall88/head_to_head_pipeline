import sys

sys.stderr = open(snakemake.log[0], "w")

import urllib.request, json

with urllib.request.urlopen(snakemake.params.mapping_url) as url:
    data = json.loads(url.read().decode())

# fix bug for one variant
bugvar = "D94H"
if bugvar in data:
    data["gyrA_" + bugvar] = data.pop(bugvar)

with open(snakemake.input.panel) as istream, open(
    snakemake.output.panel, "w"
) as ostream:
    for row in map(str.rstrip, istream):
        var_key = "_".join(row.split("\t")[:2])
        drugs = data.get(var_key)
        if not drugs:
            raise KeyError(f"No drugs found for {var_key} in mapping")
        ostream.write(row)
        print(f"\t{';'.join(drugs)}", file=ostream)
