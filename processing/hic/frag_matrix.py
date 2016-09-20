# -*- coding: utf-8 -*-
from __future__ import print_function

import argparse
import csv
import logging
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd


class Fragment(object):
    __slots__ = ["name", "start", "end", "chromosome"]

    def __init__(self, name, start, end, chromosome):
        self.name = name
        self.start = start
        self.end = end
        self.chromosome = chromosome


def load_hic_data(genome, chr1, chr2, resolution, correction, source_folder):
    file_basename = "{cell_line}_{chr1}_{chr2}_{resolution}_{correction}.txt".format(
        cell_line=genome, chr1=chr1, chr2=chr2, resolution=resolution, correction=correction)
    file_name = os.path.join(source_folder, file_basename)
    logger.info("Loading data from {file_name} into memory".format(file_name=os.path.basename(file_name)))
    result = defaultdict(dict)
    reader = csv.reader(file_name, delimiter='\t', quotechar='|')
    for row in reader:
        row_id, column_id, value = int(row[0]), int(row[1]), float(row[2])
        if row_id < column_id:
            result[row_id][column_id] = value
        else:
            result[column_id][row_id] = value
    return result


def get_fragments(fragments_file):
    logger.info("Loading fragments data from {file_name}".format(file_name=os.path.basename(fragments_file)))
    result = defaultdict(list)
    reader = csv.reader(fragments_file, delimiter="\t")
    for row in reader:
        name, start, end, chromosome = row[0], int(row[1]), int(row[2]), row[3]
        result[chromosome].append(Fragment(name=name, start=start, end=end, chromosome=chromosome))
    total = sum(lambda fr_list: len(fr_list), (result[chromosome] for chromosome in result))
    logger.info("Obtained a total of {fr_cnt} fragments over {chr_cnt} chromosomes".format(fr_cnt=total, chr_cnt=len(result)))
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
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M:%S')
    logger = logging.getLogger("frag_matrix")
    logger.setLevel(logging.DEBUG)
    fragments = get_fragments(fragments_file=args.fragments)
    contacts = pd.DataFrame(np.nan, index=fragments.keys(), columns=fragments.keys())
    for fragment1 in sorted(fragments.keys()):
        for fragment2 in sorted(fragments.keys()):
            fragments_contact = count_contact(hic_data=data, fragment1=fragments[fragment1], fragment2=fragments[fragment2], measure=args.measure)
            contacts.set_value(index=fragment1, col=fragment2, value=fragments_contact)
            contacts.set_value(index=fragment2, col=fragment1, value=fragments_contact)
    contacts.to_csv(path_or_buf=args.output)
