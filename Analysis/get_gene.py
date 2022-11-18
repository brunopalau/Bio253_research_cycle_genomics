#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 22:30:18 2022

@author: bp
"""


import pandas as pd
import numpy as np
import json

import io
import os


def read_gff(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("#")]
        
        return pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t',header=None
    ).rename(columns={0: 'seqname',1: 'source', 2: 'feature', 3:'start', 4: 'end', 5:
           'score', 6: 'strand', 7: 'frame', 8: "attribute"})
                      
                      
def get_gene(row,feature,gff):
    for line in gff.attribute[1]:
        
    
    return gff.loc[row,:]

def get_rows_from_file(path):
    output_rows = []
    cluster_rows = []
    last = False
    with open(path, "r") as f:
        for line in f.readlines():
            if last:
                output_rows.append(list(map(int,line.split(", "))))
                break
            elif line.startswith("#Dif"):
                if len(cluster_rows) != 0:
                    output_rows.append(cluster_rows)
                last = True
                continue
            elif line.startswith("#"):
                if len(cluster_rows) != 0:
                    output_rows.append(cluster_rows)
                cluster_rows = []
                continue
            elif line == "\n":
                continue
            else:
                cluster_rows.append(int(line.split(":")[0]))
    return output_rows
                
                
        
                      
# import gff
gff = read_gff("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff")
# get row we want gene from
gene_rows = get_rows_from_file("/Users/bp/Uni/Computational/HS22/BIO253/Bio253_research_cycle_genomics/Out/mutated_genes_SA6850_phenotypes.txt")
# get genes corresponding to the extracted rows
get_gene(gene_rows[1], feature, gff)
