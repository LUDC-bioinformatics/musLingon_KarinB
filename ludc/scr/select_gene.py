#!/usr/bin/python3

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "DEseq results table")
parser.add_argument("-g", "--gene", help = "potential interesting gene")
parser.add_argument("-o", "--output", help = "output")
args = parser.parse_args()

input_tab = args.input
interesting_gene = args.gene
output = open(args.output, 'w')

with open(interesting_gene, 'r') as f1:
    interesting_gene_list = []
    for line in f1:
        line = line.rstrip()
        gene = line.upper()
        interesting_gene_list.append(gene)

header = open(input_tab).readline().rstrip()
output.write(header + '\n')
with open(input_tab, 'r') as f2:
    next(f2)
    for line in f2:
        line = line.rstrip()
        GeneSymbol = line.split('\t')[-3].upper().replace('"','')
        if GeneSymbol in interesting_gene_list:
            output.write(line + '\n')
