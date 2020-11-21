#!/Users/lishuyi/miniconda3/bin/python3

import sys
import argparse
from goatools import obo_parser

parser = argparse.ArgumentParser(description="sample usage: go_parser.py -i deseq.csv -a go_anno.csv -g(optional) interested_gene_list.txt -o output.csv")
parser.add_argument("-i", "--input", help = "input DE gene list")
parser.add_argument("-a", "--annotation", help = "GO annotation list")
parser.add_argument("-g", "--gene", help = "interested gene list for searching")
parser.add_argument("-o", "--output", help = "output file name")
args = parser.parse_args()

potential_gene_list = args.input
go_annotation_list = args.annotation
interested_gene = args.gene
output = open(args.output, 'w')

term = obo_parser.GODag('reference/go-basic.obo')

with open(potential_gene_list, 'r') as f1:
    next(f1)
    gene_list = []
    gene_info = {}
    for line in f1:
        line = line.rstrip()
        gene = line.split(' ')[26].replace('"','')
        gene_list.append(gene)
        gene_info[gene] = line.split(' ')

go_dict = {}
for each in gene_list:
    go_dict[each] = []

with open(go_annotation_list, 'r') as f2:
    next(f2)
    for line in f2:
        line = line.rstrip()
        element_list = line.split(',')
        gene = element_list[0].replace('"','')
        go_dict[gene].append(element_list[1].replace('"',''))

header = open(potential_gene_list).readline().rstrip().replace(' ','\t')

output.write(header + '\t' + 'Go_id' + '\t'+ 'Go_term'+ '\t' + 'Go_num' + '\n')
for gene in go_dict:
    go_term_list = []
    for go_id in go_dict[gene]:
        if go_id in term:
            go_term_list.append(term[go_id].name)
        else:
            go_term_list.append('unidentified')
    if args.gene:
        interested_gene_list = open(interested_gene).readline().rstrip().split(',')
        if gene_info[gene][6].replace('"','') in interested_gene_list:
            output.write('\t'.join(gene_info[gene]) + '\t' + '|'.join(go_dict[gene]) + '\t' + '|'.join(go_term_list) + '\t' + str(len(go_term_list)) + '\n')
    else:
        output.write('\t'.join(gene_info[gene]) + '\t' + '|'.join(go_dict[gene]) + '\t' + '|'.join(go_term_list) + '\t' + str(len(go_term_list)) + '\n')
