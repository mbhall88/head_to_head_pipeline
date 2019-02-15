#!/usr/bin/env bash
script="$1"
out_dir="$2"
output="$3"
metadata="$4"

#download martins script to build references and then run script
perl "$script" "$out_dir"

# download lambda phage reference, make header shorter and fold lines to 60 characters
# then append it to the remove contam database
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/escherichia_virus_lambda_uid14204/NC_001416.fna -O - | \
    sed '/^>/ s/ .*//' | \
    python -c "exec(\"import sys,os\\ns=str()\\nfor l in sys.stdin.readlines():\\n\\tif l.startswith('>'): s+='{}{}{}'.format('>',l.split('|')[-2],os.linesep)\\n\\telse: s+=l.rstrip()\\nprint(s)\")" | \
    fold -w60 | \
    gzip -c -9 >> "$output"
echo -e 'Virus\t1\tNC_001416.1' >> "$metadata"
