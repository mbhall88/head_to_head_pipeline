import sys

sys.stderr = open(snakemake.log[0], "w")

import xml.etree.ElementTree as ET

tree = ET.parse(snakemake.input.xml)
root = tree.getroot()

if root.attrib["success"] != "true":
    raise ValueError("Failed response XML given")

data = []

for child in root:
    if child.tag != "SAMPLE":
        continue
    sample_acc = child.attrib["accession"]
    alias = child.attrib["alias"]
    biosample_acc = None
    for subchild in child:
        if subchild.tag != "EXT_ID":
            continue

        if biosample_acc is not None:
            raise KeyError(f"Got multiple EXT_ID tags for {alias}")
        biosample_acc = subchild.attrib["accession"]

    if biosample_acc is None:
        raise KeyError(f"Got no multiple EXT_ID tag for {alias}")

    data.append((alias, sample_acc, biosample_acc))

with open(snakemake.output.csv, "w") as fp:
    print("alias,sample,biosample", file=fp)
    for row in data:
        print(",".join(row), file=fp)
