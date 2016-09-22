# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import logging
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--frag-matrix-path", type=str, default=os.path.join("group", "cbi", "maxal", "hicproject", "software", "HiC", "processing", "hic", "frag_matrix.py"))
    parser.add_argument("--chromosomes", choices=["inter", "intra", "all"], default="all")
    parser.add_argument("--output-sh-file-prefix", default="")
    parser.add_argument("--fragments", type=str, required=True)
    parser.add_argument("--hic-dir", type=str, required=True)
    parser.add_argument("--frag-lengths", type=str, default="10000")
    parser.add_argument("--measure", choices=["inner", "outer", "fractions"], default="inner")
    parser.add_argument("--logging", type=int, choices=[logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL], default=logging.INFO)
    parser.add_argument("-o", "--output-dir", required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M:%S')
    logger = logging.getLogger("create_sh_frag_matrix_job_files")
    logger.setLevel(args.logging)
    args = parser.parse_args()
    if not os.path.exists(args.frag_matrix_path):
        logging.critical("Path {path} for frag_matrix executable does not exist".format(path=args.frag_matrix_path))
        exit(1)
    if not os.path.exists(args.fragments) or not os.path.isfile(args.fragments):
        logging.critical("Path {path} for fragments does not exist".format(path=args.fragments))
        exit(1)
    if not os.path.exists(args.hic_dir) or not os.path.isdir(args.hic_dir):
        logging.critical("Path {path} for hic files directory does not exist".format(path=args.hic_dir))
        exit(1)
    if not os.path.exists(args.output_dir):
        logger.debug("Output directory {path} does not exist. Creating one".format(path=args.output_dir))
        os.makedirs(args.output_dir)
    file_template = "\n".join([
        "#!/bin/sh",
        "#SBATCH -p short",
        "#SBATCH -J {base_name}",
        "#SBATCH -t 2:00:00",
        "#SBATCH -e {base_name}.err",
        "#SBATCH -o {base_name}.out",
        "module load python/2.7.6",
        "python {frag_matrix_path} --fragments {fragments_path} --hic {hic_path} --frag-lengths {frag_lengths} --measure {measure} --output {output_file_rel_path}"
    ])
    hic_export_files = [f for f in os.listdir(args.hic_dir) if f.endswith(".txt") and len(f.split("_")) > 4]
    logger.info("Found {hic_cnt} hic files".format(hic_cnt=len(hic_export_files)))
    for hic_file in hic_export_files:
        logger.debug("Creating batch sh file for {hic_file_name}".format(hic_file_name=hic_file))
        cell_line, chr1, chr2, resolution, correction = hic_file.split("_")
        if chr1 == chr2 and args.chromosomes == "inter":
            continue
        if chr1 != chr2 and args.chromosomes == "intra":
            continue
        if correction.endswith(".txt"):
            correction = correction[:-4]
        sh_file_name = "{prefix}{chr1}_{chr2}.sh".format(chr1=chr1, chr2=chr2, prefix=args.output_sh_file_prefix)
        contacts_file_name = "all_{frag_lengths}_{chr1}_{chr2}_{bin_size}_KR.txt".format(
            frag_lengths=args.frag_lengths, chr1=chr1, chr2=chr2, bin_size=args.bin_size)
        with open(os.path.join(args.output_dir, sh_file_name), "wt") as dest:
            print(file_template.format(base_name=hic_file[:-4],
                                       frag_matrix_path=os.path.abspath(args.frag_matrix_path),
                                       fragments_path=os.path.abspath(args.fragments),
                                       hic_path=os.path.join(args.hic_dir, hic_file),
                                       frag_lengths=args.frag_lengths,
                                       measure=args.measure,
                                       output_file_rel_path=contacts_file_name), file=dest)
