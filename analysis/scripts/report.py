import json
import os
import sys
from snakemake.utils import report


# capture all stderr messages in log file
class Logger(object):
    """Class to capture stderr in log file"""
    def __init__(self, log_path):
        self.log = open(log_path, "w")

    def write(self, message):
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

sys.stderr = Logger(snakemake.log[0])


def mykrobe_overview(filepath: str) -> dict:
    """Extracts the susceptiblity information from mykrobe predict json.
    :param filepath: path to mykrobe predict output json file.
    :returns A dictionary of the susceptibility results.
    """
    sample_id = os.path.basename(filepath).split('.')[0].split('_predict')[0]
    with open(filepath, 'r') as mykrobe_json:
        data = json.load(mykrobe_json)
    return data[sample_id].get('susceptibility', {})

def mykrobe_rst_list(data: dict) -> str:
    """Formats the Mykrobe data into a restructuredtext unordered list."""
    result = ''
    for drug in data:
        drug_info = data.get(drug)
        predict = drug_info.get('predict', '')
        if predict == 'S':
            result += '- **{drug}**\n    - Prediction: **Susceptible**\n'.format(drug=drug)
        elif predict == 'R':
            called_by = list(drug_info.get('called_by', []))
            result += '- **{drug}**\n    - Prediction: **Resistant**\n'.format(drug=drug)
            for var in called_by:
                coverage = drug_info.get('called_by').get(var).get('info', '').get('coverage', '')
                ref_coverage = coverage.get('reference', '').get('median_depth', '')
                alt_coverage = coverage.get('alternate', '').get('median_depth', '')
                result += '    - Called by: {var}\n'.format(var=var)
                result += '    - Reference median depth: {ref_coverage}\n'.format(ref_coverage=ref_coverage)
                result += '    - Alternate median depth: {alt_coverage}\n'.format(alt_coverage=alt_coverage)
    return result

def get_phylo_group_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("phylo_group", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_species_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("species", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_lineage_string(d):
    s = []
    depth=[]
    per_cov=[]
    for k, v in d.get("phylogenetics", {}).get("lineage", {}).items():
        s.append(k)
        depth.append(str(v.get("median_depth")))
        per_cov.append(str(v.get("percent_coverage")))
    return ";".join(s), ";".join(per_cov), ";".join(depth)


def get_num_reads(stats_file: str) -> int:
    """Extracts the number of reads from a nanostats text file."""
    with open(stats_file, 'r') as stats:
        for line in stats:
            if 'Number of reads:' in line:
                num_reads = line.split()[-1]
    return int(num_reads)


def load_json(filepath: str) -> dict:
    with open(filepath, 'r') as infile:
        data = json.load(infile)
    return data


mykrobe_report = mykrobe_rst_list(mykrobe_overview(snakemake.input.mykrobe))
num_reads_pre_filter = get_num_reads(snakemake.input.stats_pre_filter)
num_reads_post_filter = get_num_reads(snakemake.input.stats_post_filter)
percent_reads_mapped = round(num_reads_post_filter / num_reads_pre_filter * 100, 2)
sample=snakemake.params.sample

mykrobe_data = load_json(snakemake.input.mykrobe)
key = list(mykrobe_data.keys())[0]
data = mykrobe_data[key]
phylo_group, _, __ = get_phylo_group_string(data)
species, _, __ = get_species_string(data)
lineage, _, __ = get_lineage_string(data)

with open(snakemake.input.trim_log, 'r') as log_file:
    trim_log = log_file.read()


report("""
===================================
Report for {sample}
===================================
Processing Steps
===================================
1. Fast5 files were basecalled with Guppy.
2. Resulting fastq files were demultiplexed with Deepbinner_ (default settings). More detailed information about the demultiplexing can be found in `demultiplex_log`_. For statistics of each sample at this stage, see `stats_pre_filter`_.
3. Porechop_ was run to trim adapter sequences from reads and discard reads with adapters found in the middle. More detailed information can be found in `trim_log`_. For quality control plots of the reads after this step, see `qc_plot`_.
4. Reads were aligned to the CRyPTIC TB contamination database using Minimap2_.
5. Plotting of the sample compostion after alignment to this database was performed with Krona_. The plot can be found in `compostion_plot`_.
6. All reads which did not map to TB or NTMs were filtered out. Prior to filtering there were {num_reads_pre_filter} reads. After filtering there remains {num_reads_post_filter}. This means {percent_reads_mapped}% of reads mapped to TB/NTMs. For more statistics on the post-filtered reads see `stats_post_filter`_. For quality control plots of the reads after this step (and read percent identity to the database) see `qc_plot`_. Stats were produced with NanoStat_ and plots with Pistis_.
7. `Mykrobe predict`_ was run on the filtered fastq files to predict drug susceptiblity.

Mykrobe Analysis
===================================
**Phylogenetic group:** {phylo_group}
**Species:** {species}
**Lineage:** {lineage}
A summary of the susceptiblity information from `Mykrobe predict`_ is shown here. For the full report, see mykrobe_. If resistance is identified for a drug then the predicted responsible variant(s) is given, along with supporting information.
{mykrobe_report}
.. _Porechop: https://github.com/rrwick/Porechop
.. _Minimap2: https://github.com/lh3/minimap2
.. _NanoStat: https://github.com/wdecoster/nanostat
.. _Pistis: https://github.com/mbhall88/pistis
.. _`Mykrobe predict`: https://github.com/Mykrobe-tools/mykrobe
.. _Deepbinner: https://github.com/rrwick/Deepbinner
.. _Krona: https://github.com/marbl/Krona
""", snakemake.output[0], metadata="Author: Michael Hall (michael.hall@ebi.ac.uk)",
**snakemake.input)
