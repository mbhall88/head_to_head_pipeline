#!/usr/bin/env python3
"""Pairwise SNP distance bwtween all pair combinations between two fasta files.
The files must have the same number of sequences, the same lengths for all sequences
and the same fasta head IDs. The output is a CSV distance matrix where the rows
represent <a.fa> and columns <b.fa>"""
import sys
from itertools import product, starmap
from pathlib import Path
from typing import Collection

import numpy as np
import pyfastx
from tqdm import tqdm

N = "N"
DELIM = ","
CLEAN_IDX = True


def is_same(x: str, y: str) -> bool:
    return N in (x, y) or x == y


def is_diff(x: str, y: str) -> bool:
    return not is_same(x, y)


def hamming(s1: Collection, s2: Collection) -> int:
    assert len(s1) == len(s2)
    return int(sum(starmap(is_diff, zip(s1, s2))))


def main():
    if len(sys.argv) < 3:
        print(f"USAGE: {sys.argv[0]} <a.fa> <b.fa>")
        sys.exit(1)

    fp1 = sys.argv[1]
    fp2 = sys.argv[2]
    if CLEAN_IDX:
        for fp in (fp1, fp2):
            Path(fp + ".fxi").unlink(missing_ok=True)

    fa1 = pyfastx.Fasta(fp1)
    fa2 = pyfastx.Fasta(fp2)

    n_seqs = len(fa1)

    if n_seqs != len(fa2):
        raise ValueError("Different number of sequences in the two files")

    ids = sorted(fa1.keys())

    all_seq_ids_same = all(x == y for x, y in zip(ids, sorted(fa2.keys())))
    if not all_seq_ids_same:
        raise ValueError("Sequence IDs are not identical between the two files")

    mtx = np.zeros(shape=(n_seqs, n_seqs), dtype=np.int)

    for i, j in tqdm(product(range(n_seqs), range(n_seqs)), total=n_seqs * n_seqs):
        dist = hamming(fa1[ids[i]].seq, fa2[ids[j]].seq)
        mtx[i][j] = dist

    print(DELIM.join(["sample", *ids]))
    for i, sample in enumerate(ids):
        row = DELIM.join(map(str, mtx[i]))
        print(f"{sample}{DELIM}{row}")

    if CLEAN_IDX:
        for fp in (fp1, fp2):
            Path(fp + ".fxi").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
