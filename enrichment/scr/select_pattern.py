#!/Users/lishuyi/miniconda3/bin/python3

import sys
import argparse

parser = argparse.ArgumentParser(description="sample usage: select_pattern.py -i GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab -p Pattern1 -o pattern1.csv")
parser.add_argument("-i", "--input", help = "RSEM results table")
parser.add_argument("-p", "--pattern", help = "pattern")
parser.add_argument("-o", "--output", help = "output")
args = parser.parse_args()

input_tab = args.input
filter = args.pattern
output = open(args.output, 'w')

output.write('EnsemblID' + '\t' + 'GeneSymbol' + '\n')
with open(input_tab, 'r') as f1:
    next(f1)
    for line in f1:
        line = line.rstrip()
        pattern = line.split('\t')[6].replace('"','')
        gene = line.split('\t')[0].replace('"','')
        EnsemblID = gene.split('_')[0]
        GeneSymbol = gene.split('_')[1]
        if pattern == filter:
            output.write(EnsemblID + '\t' + GeneSymbol + '\n')
