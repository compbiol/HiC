#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import sys
import networkx as nx
import csv
from six import next

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mapping", nargs="+", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("--skip_header", action="store_true", default=False)
    parser.add_argument("--separator", default="\t")
    args = parser.parse_args()
    graph = nx.Graph()
    for input_file in args.mapping:
        reader = csv.reader(input_file, delimiter=args.separator)
        next(reader)
        for row in reader:
            graph.add_edge(row[0], row[1])

    for cnt, cc in enumerate(nx.connected_component_subgraphs(graph)):
        for node in cc:
            print(node, cnt, file=args.output)
