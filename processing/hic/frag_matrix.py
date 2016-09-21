#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

import argparse
import csv
import datetime
import logging
import os
import sys
import math
from collections import defaultdict


class Fragment(object):
    __slots__ = ["name", "start", "end", "chromosome"]

    def __init__(self, name, start, end, chromosome):
        self.name = name
        self.start = start
        self.end = end
        self.chromosome = chromosome


def load_hic_data(hic_filename):
    result = defaultdict(dict)
    with open(hic_filename, "rt") as source:
        reader = csv.reader(source, delimiter='\t', quotechar='|')
        for cnt, row in enumerate(reader):
            if row[0].strip().startswith("#"):
                logger.debug("Skipping row {row}".format(row="\t".join(row)))
                continue
            row_id, column_id, value = int(row[0]), int(row[1]), float(row[2])
            row, column = (row_id, column_id) if row_id < column_id else (column_id, row_id)
            result[row][column] = value
    return result


def get_fragments(fragments_filename):
    result = []
    with open(fragments_filename, "rt") as source:
        reader = csv.reader(source, delimiter="\t")
        for row in reader:
            name, start, end, chromosome = row[0], int(row[1]), int(row[2]), row[3]
            result.append(Fragment(name=name, start=start, end=end, chromosome=chromosome))
    return result


def count_contacts(hic_data, f1, f2, measure, step):
    def get_contact(data, i, j, im, jm, observed):
        index1, index2 = (i, j) if i < j else (j, i)
        if (index1, index2) in observed:
            return 0
        contact = data[index1].get(index2, 0) * im * jm
        observed.add((index1, index2))
        return contact

    def get_before_after_indexes(f, fsa, feb, fully_inside, measure):
        if measure == "outer":
            before, after = [1], [1]
        elif measure == "inner":
            before, after = [0], [0]
        else:  # fractions
            before, after = [(fsa - f.start) / step], [(f.end - feb) / step]
        if fully_inside:
            after = []
            before = [((f.start - feb) + (fsa - f.end)) / step], []
        return before, after

    result = 0
    f1sb = int(math.floor(f1.start / step) * step)
    f1sa = int(math.ceil(f1.start / step) * step)
    f1eb = int(math.floor(f1.end / step) * step)
    f1ea = int(math.ceil(f1.end / step) * step)

    if f1eb <= f1sa:
        logger.warning("Fragment {fragment} is short ({f_length} is less than 2 bins of size {size}).".format(fragment=f1.name, size=step,
                                                                                                              f_length=f1.end - f1.start))
        logger.debug("Fragment {f}: start={start}, end={end}, fsb={fsb}, fsa={fsa}, feb={feb}, fea={fea}"
                     "".format(f=f1.name, start=f1.start, end=f1.end, fsb=f1sb, fsa=f1sa, feb=f1eb, fea=f1ea))

    f2sb = int(math.floor(f2.start / step) * step)
    f2sa = int(math.ceil(f2.start / step) * step)
    f2eb = int(math.floor(f2.end / step) * step)
    f2ea = int(math.ceil(f2.end / step) * step)

    if f2eb <= f2sa:
        logger.warning("Fragment {fragment} is short ({f_length} is less than 2 bins of size {size}).".format(fragment=f2.name, size=step,
                                                                                                              f_length=f2.end - f2.start))
        logger.debug("Fragment {f}: start={start}, end={end}, fsb={fsb}, fsa={fsa}, feb={feb}, fea={fea}"
                     "".format(f=f2.name, start=f2.start, end=f2.end, fsb=f2sb, fsa=f2sa, feb=f2eb, fea=f2ea))

    f1_inner_multipliers = [1] * (f1eb - f1sa)  # in case of short fragment multiplication will happen with negative number and en empty list will be produced
    f2_inner_multipliers = [1] * (f2eb - f2sa)  # same

    f1b, f1a = get_before_after_indexes(f=f1, fsa=f1sa, feb=f1eb, fully_inside=f1sb == f1eb and f1sa == f1ea, measure=measure)
    f2b, f2a = get_before_after_indexes(f=f1, fsa=f1sa, feb=f1eb, fully_inside=f2sb == f2eb and f2sa == f2ea, measure=measure)

    f1_multipliers = f1b + f1_inner_multipliers + f1a
    f2_multipliers = f2b + f2_inner_multipliers + f2a

    f1_indexes = range(f1sb, f1ea, step)
    f2_indexes = range(f2sb, f2ea, step)

    counted = set()
    for f1_index, f1_multiplier in zip(f1_indexes, f1_multipliers):
        for f2_index, f2_multiplier in zip(f2_indexes, f2_multipliers):
            result += get_contact(data=hic_data, i=f1_index, im=f1_multiplier, j=f2_index, jm=f2_multiplier, observed=counted)
    return result


def get_chromosomes_from_hic_filename(hic):
    # file name template "cell-line_chr1_chr2_resolution_correction.txt"
    return os.path.basename(hic).split("_")[1:3]


def get_step_from_hic_filename(hic):
    return int(os.path.basename(hic).split("_")[3])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hic", type=str, required=True)
    parser.add_argument("--fragments", type=str, required=True)
    parser.add_argument("--measure", choices=["outer", "inner", "fractions"], default="inner")
    parser.add_argument("--existing", type=str, default=None)
    parser.add_argument("--frag-lengths", type=int, default=-1)
    parser.add_argument("--logging", type=int, choices=[logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL], default=logging.INFO)
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M:%S')
    logger = logging.getLogger("frag_matrix")
    logger.setLevel(args.logging)
    start_time = datetime.datetime.now()
    logger.info("Starting the whole show...")
    logger.info("Loading fragments info from {f_filename}".format(f_filename=os.path.basename(args.fragments)))
    fragments = get_fragments(fragments_filename=args.fragments)
    logger.info("Loaded a total of {f_cnt} fragments".format(f_cnt=len(fragments)))
    if args.frag_lengths > 0:
        logger.info("Filtering out fragments shorter than {f_length_thresh}".format(f_length_thresh=args.frag_lengths))
        fragments = [f for f in fragments if f.end - f.start > args.frag_lengths]
        logger.info("A total of {f_cnt} are longer than {f_length_thresh}".format(f_cnt=len(fragments), f_length_thresh=args.frag_lengths))
    chr1, chr2 = get_chromosomes_from_hic_filename(hic=args.hic)
    step = get_step_from_hic_filename(hic=args.hic)
    logger.info("Working with chromosomes {chr1} and {chr2} and a step of {step}".format(chr1=chr1, chr2=chr2, step=step))
    logger.info("Filtering fragments that do not belong to {chr1} or {chr2} chromosomes".format(chr1=chr1, chr2=chr2))
    fragments = [fragment for fragment in fragments if fragment.chromosome in [chr1, chr2]]
    if chr1 == chr2:
        fragments1 = fragments[:]
        fragments2 = fragments[:]
    else:
        fragments1 = [f for f in fragments if f.chromosome == chr1]
        fragments2 = [f for f in fragments if f.chromosome == chr2]
    logger.info("Will be computing pairwise contact between {f1_cnt} and {f2_cnt} fragments".format(f1_cnt=len(fragments1), f2_cnt=len(fragments2)))
    logger.info("Loading HiC data from {hic_filename}".format(hic_filename=os.path.basename(args.hic)))
    data = load_hic_data(hic_filename=args.hic)
    contacts = defaultdict(dict)
    observed = set()
    logger.info("Computing contacts")
    for fragment1 in fragments1:
        for fragment2 in fragments2:
            f1, f2 = (fragment1.name, fragment2.name) if fragment1.name < fragment2.name else (fragment2.name, fragment1.name)
            if (f1, f2) in observed:
                continue
            logger.debug("Computing contacts between fragments {f1} and {f2}".format(f1=fragment1.name, f2=fragment2.name))
            fragments_contact = count_contacts(hic_data=data, f1=fragment1, f2=fragment2, measure=args.measure, step=step)
            logger.debug("Done. {f1} and {f2} have {c_cnt} HiC contacts".format(f1=fragment1.name, f2=fragment2.name, c_cnt=fragments_contact))
            observed.add((f1, f2))
            contacts[f1][f2] = fragments_contact
    logger.info("Computed all pairwise contacts. Outputting results.")
    print("# python :: {python_version}".format(python_version=".".join(map(str, sys.version_info))))
    print("# hic source :: {hic} ".format(hic=os.path.abspath(args.hic)), file=args.output)
    print("# fragments source :: {fragments}".format(fragments=os.path.abspath(args.fragments)), file=args.output)
    print("# measure :: {measure}".format(measure=args.measure), file=args.output)
    print("# bin size :: {step}".format(step=step), file=args.output)
    print("# fragment min. size :: {f_lengths}".format(f_lengths=args.frag_lengths), file=args.output)
    print("# chromosome 1 :: {f1}".format(f1=chr1), file=args.output)
    print("# chromosome 2 :: {f2}".format(f2=chr2), file=args.output)
    print("# chromosome {chr1} fragment cnt :: {f1_cnt}".format(f1_cnt=len(fragments1), chr1=chr1), file=args.output)
    print("# chromosome {chr2} fragment cnt :: {f2_cnt}".format(f2_cnt=len(fragments2), chr2=chr2), file=args.output)
    for key1 in sorted(contacts.keys()):
        for key2 in sorted(contacts[key1].keys()):
            print(key1, key2, contacts[key1][key2], sep="\t", file=args.output)
    logger.info("All done. It only took us: {time_cnt}".format(time_cnt=str(datetime.datetime.now() - start_time)))
