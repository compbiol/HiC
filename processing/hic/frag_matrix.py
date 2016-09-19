# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import sys
from collections import defaultdict

import pandas as pd
import numpy as np
import os


def get_hic_data(genome, source_folder):
    if source_folder is None:
        source_folder = genome
    file_names = [file_name for file_name in os.listdir(genome) if file_name.startswith(genome)]

    result = defaultdict(dict)
    return result


def get_fragments(fragments_file):
    result = {}
    return result


def count_contact(hic_data, fragment1, fragment2, measure):
    result = np.nan
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=str, required=True)
    parser.add_argument("--hic-source", type=str, default=None)
    parser.add_argument("--fragments", type=argparse.FileType("rt"), required=True)
    parser.add_argument("--measure", choices=["outer", "inner", "fraction"], default="outer")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()
    data = get_hic_data(genome=args.genome, source_folder=args.hic_source)
    fragments = get_fragments(fragments_file=args.fragments)
    contacts = pd.DataFrame(np.nan, index=fragments.keys(), columns=fragments.keys())
    for fragment1 in sorted(fragments.keys()):
        for fragment2 in sorted(fragments.keys()):
            fragments_contact = count_contact(hic_data=data, fragment1=fragments[fragment1], fragment2=fragments[fragment2], measure=args.measure)
            contacts.set_value(index=fragment1, col=fragment2, value=fragments_contact)
            contacts.set_value(index=fragment2, col=fragment1, value=fragments_contact)
    contacts.to_csv(path_or_buf=args.output)
