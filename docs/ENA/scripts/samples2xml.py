"""
This scripts takes a samplesheet and creates an ENA-compliant Sample XML file.
See https://ena-docs.readthedocs.io/en/latest/submit/samples/programmatic.html and
https://ena-docs.readthedocs.io/en/latest/submit/samples.html
USAGE: python samples2xml.py samplesheet.csv
"""
# Good example for Sample entry https://www.ebi.ac.uk/ena/browser/view/SAMEA95444668
import sys
import xml.etree.ElementTree as ET

import pandas as pd

# https://www.ebi.ac.uk/ena/browser/view/ERC000028
CHECKLIST = "ERC000028"
TAXON_ID = ET.Element("TAXON_ID")
TAXON_ID.text = "1773"
SCI_NAME = ET.Element("SCIENTIFIC_NAME")
SCI_NAME.text = "Mycobacterium tuberculosis"
TITLE = ET.Element("TITLE")
TITLE.text = "Mycobacterium tuberculosis WGS"
LOCATION = {
    "madagascar": "Madagascar",
    "south_africa": "South Africa",
    "birmingham": "United Kingdom",
}


def attribute(tag: str, value) -> ET.Element:
    el = ET.Element("SAMPLE_ATTRIBUTE")
    tag_el = ET.Element("TAG")
    tag_el.text = tag
    el.append(tag_el)
    value_el = ET.Element("VALUE")
    value_el.text = value
    el.append(value_el)
    return el


def attributes(date: str, loc: str, lineage: str) -> ET.Element:
    attrs = ET.Element("SAMPLE_ATTRIBUTES")
    attrs.append(attribute("collection date", date))
    attrs.append(attribute("strain", f"Lineage {lineage}"))
    attrs.append(attribute("host health state", "diseased"))
    attrs.append(attribute("host scientific name", "homo sapien"))
    attrs.append(attribute("isolate", "Mycobacterium tuberculosis"))
    attrs.append(attribute("isolation_source", "sputum"))
    attrs.append(attribute("geographic location (country and/or sea)", loc))
    attrs.append(attribute("ENA-CHECKLIST", CHECKLIST))

    return attrs


def add_sample_name(el: ET.Element):
    name = ET.Element("SAMPLE_NAME")
    name.append(TAXON_ID)
    name.append(SCI_NAME)
    el.append(name)


def main():
    if len(sys.argv) < 2:
        raise ValueError("No samplesheet passed")
    samplesheet = pd.read_csv(sys.argv[1])
    print('<?xml version="1.0" encoding="UTF-8"?>')
    root = ET.Element("SAMPLE_SET")
    for _, row in samplesheet.iterrows():
        sample_el = ET.Element("SAMPLE", attrib=dict(alias=row["sample"]))
        sample_el.append(TITLE)
        add_sample_name(sample_el)
        location = LOCATION[row.site]
        date = row["collection_date"]
        if pd.isna(date):
            print(
                f"[WARNING] Date is null for sample {row['sample']}...skipping...",
                file=sys.stderr,
            )
            continue
        attrs = attributes(date, location, row.lineage)
        sample_el.append(attrs)
        root.append(sample_el)

    ET.indent(root)
    ET.dump(root)


if __name__ == "__main__":
    main()
