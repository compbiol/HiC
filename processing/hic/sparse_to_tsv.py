#! /usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import csv
import logging
import sys
import pandas as pd
import numpy as np


def read_hic_export(file_name, step=None, ignore_empty_start=True, symmetrical=True):
    data = []
    chr1_max_value = -1
    chr2_max_value = -1
    chr1_min_value = 100000000000
    chr2_min_value = 100000000000
    reader = csv.reader(file_name, delimiter='\t', quotechar='|')
    for row in reader:
        row_id, column_id, value = int(row[0]), int(row[1]), float(row[2])
        data.append((row_id, column_id, value))
        if column_id > chr2_max_value:
            chr2_max_value = column_id
        if row_id > chr1_max_value:
            chr1_max_value = row_id
        if column_id < chr2_min_value:
            chr2_min_value = column_id
        if row_id < chr1_min_value:
            chr1_min_value = row_id
    if step is None:
        rows = sorted(set([entry[0] for entry in data]))
        rows_min = min(a2 - a1 for a1, a2 in zip(rows[:-1], rows[1:]))
        columns = sorted(set([entry[1] for entry in data]))
        columns_min = min(a2 - a1 for a1, a2 in zip(columns[:-1], columns[1:]))
        step = min(rows_min, columns_min)
        logger.info("Inferred matrix resolution: {res}".format(res=step))
    chr1_min_value = 0 if ignore_empty_start else chr1_min_value
    chr2_min_value = 0 if ignore_empty_start else chr2_min_value
    if symmetrical:
        if chr1_min_value < chr2_min_value:
            chr2_min_value = chr1_min_value
        else:
            chr1_min_value = chr2_min_value
        if chr1_max_value > chr2_max_value:
            chr2_max_value = chr1_max_value
        else:
            chr1_max_value = chr2_max_value
    logger.info("Chr1 (rows) range from: {start} to {end} with step {step}".format(start=chr1_min_value,
                                                                                   end=chr1_max_value,
                                                                                   step=step))
    logger.info("Chr2 (columns) range from: {start} to {end} with step {step}".format(start=chr2_min_value,
                                                                                      end=chr2_max_value,
                                                                                      step=step))
    df = pd.DataFrame(0, index=np.arange(chr1_min_value, chr1_max_value, step),
                      columns=np.arange(chr2_min_value, chr2_max_value, step))
    for row_id, column_id, value in data:
        value = 0.0 if np.isnan(value) else value
        df.set_value(row_id, column_id, value)
        df.set_value(column_id, row_id, value)

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sparse_contact_matrix", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("--ignore-empty-start", action="store_true", default=False)
    parser.add_argument("--step", default=None)
    parser.add_argument("--diff-chromosomes", action="store_true", dest="same_chromosomes", default=False)
    parser.add_argument("--output-separator", default="\t")
    parser.add_argument("--output-upper-triangular", action="store_true", dest="upper_tria", default=False)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    logger = logging.getLogger("sparse_to_csv")
    logger.setLevel(logging.DEBUG)
    logger.info("Processing file {file_name}".format(file_name=args.sparse_contact_matrix.name))
    if args.step is None:
        logger.info("Matrix resolution is not specified, will be inferred")
    else:
        logger.info("Matrix resolution is specified at {res}".format(res=args.step))
    df = read_hic_export(file_name=args.sparse_contact_matrix,
                         step=args.step,
                         ignore_empty_start=args.ignore_empty_start,
                         symmetrical=args.same_chromosomes)
    df = df.fillna(0.0)
    if args.upper_tria:
        logger.info("Substituting all the data from lower triangle of the matrix with zeros")
        df = df.where(np.triu(np.ones(df.shape)).astype(np.bool)).fillna(0)
    df.to_csv(path_or_buf=args.output, sep=args.output_separator)
