# -*- coding: utf-8 -*-
from __future__ import print_function, division

import argparse
import logging

import sys

HUMAN_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X"]
MOUSE_CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX"]
DOG_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "X"]
RABBIT_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X"]
MACAQUE_CHROMOSOMES = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "X"]

CELL_LINES_CHROMOSOMES = {
    "rao14_gm12878": HUMAN_CHROMOSOMES,
    "rao14_imr90": HUMAN_CHROMOSOMES,
    "rao14_hmec": HUMAN_CHROMOSOMES,
    "rao14_nhek": HUMAN_CHROMOSOMES,
    "rao14_k562": HUMAN_CHROMOSOMES,
    "rao14_kbm7": HUMAN_CHROMOSOMES,
    "rao14_huvec": HUMAN_CHROMOSOMES,
    "rao14_chx12": MOUSE_CHROMOSOMES,
    "rudan15_dog1": DOG_CHROMOSOMES,
    "rudan15_dog2": DOG_CHROMOSOMES,
    "rudan15_mouse1": MOUSE_CHROMOSOMES,
    "rudan15_mouse2": MOUSE_CHROMOSOMES,
    "rudan15_macaque1": MACAQUE_CHROMOSOMES,
    "rudan15_macaque2": MACAQUE_CHROMOSOMES,
    "rudan15_rabbit1": RABBIT_CHROMOSOMES,
    "rudan15_rabbit2": RABBIT_CHROMOSOMES
}

CELL_LINE_SOURCES = {
    "rao14_gm12878": "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic",
    "rao14_imr90": "https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic",
    "rao14_hmec": "https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic",
    "rao14_nhek": "https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic",
    "rao14_k562": "https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic",
    "rao14_kbm7": "https://hicfiles.s3.amazonaws.com/hiseq/kbm7/in-situ/combined.hic",
    "rao14_huvec": "https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic",
    "rao14_hela": "https://hicfiles.s3.amazonaws.com/hiseq/hela/in-situ/combined.hic",
    "rao14_chx12": "https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic",
    "rudan15_dog1": "http://hicfiles.s3.amazonaws.com/external/rudan/canis-lupus-rep1.hic",
    "rudan15_dog2": "http://hicfiles.s3.amazonaws.com/external/rudan/canis-lupus-rep1.hic",
    "rudan15_mouse1": "http://hicfiles.s3.amazonaws.com/external/rudan/mouse-rep1.hic",
    "rudan15_mouse2": "http://hicfiles.s3.amazonaws.com/external/rudan/mouse-rep2.hic",
    "rudan15_macaque1": "http://hicfiles.s3.amazonaws.com/external/rudan/macaque-rep1.hic",
    "rudan15_macaque2": "http://hicfiles.s3.amazonaws.com/external/rudan/macaque-rep2.hic",
    "rudan15_rabbit1": "http://hicfiles.s3.amazonaws.com/external/rudan/rabbit-rep1.hic",
    "rudan15_rabbit2": "http://hicfiles.s3.amazonaws.com/external/rudan/rabbit-rep2.hic"
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell-line", type=str, choices=CELL_LINES_CHROMOSOMES.keys(), required=True)
    parser.add_argument("--chromosomes", type=str, choices=["inter", "intra", "all"], required=True)
    parser.add_argument("--data", type=str, choices=["observed", "oe", "pearson", "norm", "expected"], default="observed")
    parser.add_argument("--norm", type=str, choices=["NONE", "VC", "VC_SQRT", "KR"], default="KR")
    parser.add_argument("--bin-size", type=str, default="5000")
    parser.add_argument("--juicebox-tools-path", type=str, required=True)
    parser.add_argument("--job-name", default="", type=str)
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    logger = logging.getLogger("creating_sh_export_file")
    args = parser.parse_args()
    chromosomes = []
    if args.chromosomes == "inter" or args.chromosomes == "all":
        observed = set()
        for chr1 in CELL_LINES_CHROMOSOMES[args.cell_line]:
            for chr2 in CELL_LINES_CHROMOSOMES[args.cell_line]:
                c1, c2 = (chr1, chr2) if chr1 < chr2 else (chr2, chr1)
                if (c1, c2) in observed:
                    continue
                if c1 == c2 and args.chromosomes == "inter":
                    continue
                chromosomes.append((c1, c2))
                observed.add((c1, c2))
    else:
        chromosomes = [(c, c) for c in CELL_LINES_CHROMOSOMES[args.cell_line]]
    job_name = args.job_name
    if len(args.job_name) == 0:
        job_name = "_".join([args.cell_line.replace("_", "-"), args.chromosomes, args.data, args.norm, args.bin_size])
    print("#!/bin/sh", file=args.output)
    print("#SBATCH -p short", file=args.output)
    print("#SBATCH -t 6:00:00")
    print("#SBATCH -J {job_name}".format(job_name=job_name), file=args.output)
    print("#SBATCH -o {job_name}.out".format(job_name=job_name), file=args.output)
    print("#SBATCH -e {job_name}.err".format(job_name=job_name), file=args.output)
    print("module load jdk/1.8.0", file=args.output)
    for chr1, chr2 in chromosomes:
        c1 = chr1[3:] if chr1.startswith("chr") else chr1
        c2 = chr2[3:] if chr2.startswith("chr") else chr2
        output_path = "{cell_line}_{chr1}_{chr2}_{bin_size}_{correction}.txt".format(cell_line=args.cell_line.replace("_", "-"),
                                                                                     chr1=c1, chr2=c2, bin_size=args.bin_size,
                                                                                     correction=args.norm)
        print("echo \"working with cell-line: {cell_line}; chromosome 1: {chr1}; chromosomes 2: {chr2}; resolution {bin_size}\""
              "".format(cell_line=args.cell_line, chr1=chr1, chr2=chr2, bin_size=args.bin_size), file=args.output)
        print("java -jar {juicebox_tools_path} dump {data} {correction} {path} {chr1} {chr2} BP {bin_size} {output_path}"
              "".format(juicebox_tools_path=args.juicebox_tools_path,
                        data=args.data, correction=args.norm, path=CELL_LINE_SOURCES[args.cell_line],
                        chr1=chr1, chr2=chr2, bin_size=args.bin_size, output_path=output_path),
              file=args.output)
