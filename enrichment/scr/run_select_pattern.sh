#!/bin/bash

mkdir Matrix/pattern

for i in $(seq -w 1 5)
do
  scr/select_pattern.py -i Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults -p Pattern${i} -o Matrix/pattern/pattern${i}.csv
  scr/select_pattern.py -i Matrix/GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab -p Pattern${i} -o Matrix/pattern/pattern${i}_FDR0.05.csv
done
