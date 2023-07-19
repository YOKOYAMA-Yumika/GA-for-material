#!/usr/bin/env python
# coding: utf-8

import os
import sys

elem_list = [[38, 26], [13, 51]]

strings = []
for elem in elem_list:
    gene = ""
    for x, i in enumerate(elem):
        for _ in range(0, i):
            gene += str(x)
    strings.append(gene)

with open("inp_ga.py", "r") as f:
    inpga = f.readlines()
inpga_n = []

strings_join = " ".join(strings)
for i in inpga:
    if "GENOMS" in i:
        i = f'GENOMS = "{strings_join}" #original gen(use as reference) str'
    inpga_n.append(i)

with open("inp.params", "w") as w:
    for i in strings:
        w.write(i + "\n")
    
