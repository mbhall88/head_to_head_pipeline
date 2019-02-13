#!/usr/bin/env sh

is_multiplexed="$1"
fast5_dir="$2"
classification_path="$3"
fastq="$4"
out_dir="$5"
output="$6"

if [[ $is_multiplexed = "True" ]]; then
    deepbinner classify --native "$fast5_dir" > "$classification_path"

    echo "Deepbinner classification finished."

    deepbinner bin --classes "$classification_path" \
        --reads "$fastq" \
        --out_dir "$out_dir"
fi

touch "$output"
