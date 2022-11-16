#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:09:54 2022

@author: bp
"""

import pandas as pd
import numpy as np


# import vcf file
import io
import os

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})



def change_name(df):
    new_names = []
    for id_,cond in zip(df.ID,df.Condition):
        
        new_string = "S68_"
        
        if cond == "supernatant":
            new_string += "S_"
        else:
            new_string += "T_"
        
        letter = id_.split("_")[1][0]
        
        num1,num2 = id_.split(letter)[1].split(".")
        
        if len(num1) == 1:
            num1 = "0" + num1
        
        new_string += letter
        new_string += num1
        new_string += "_"
        new_string += num2
        
        new_names.append(new_string)
    df["new_names"] = new_names
    return df


def read_gff(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("#")]
        
        return pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t',header=None
    ).rename(columns={0: 'seqname',1: 'source', 2: 'feature', 3:'start', 4: 'end', 5:
           'score', 6: 'strand', 7: 'frame', 8: "attribute"})
    
                      
# import clustering results
df = pd.read_csv('/Users/bp/Uni/Computational/HS22/BIO253/Data/cluster_results.csv')
# adjust names
df = change_name(df)
# import vcf
vcf = read_vcf("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/2020-09-25_mod_SA6850_filteredFromAncestor.vcf")


# list of id for each cluster
cluster_1_id = df[df["phenotype_clusters_2"] == 1]
cluster_2_id = df[df["phenotype_clusters_2"] == 2]


# create dicitonaries to later safe mutations
mut_in_1 = {}
mut_in_2 = {}

def get_mut(id_list, vcf_df):
    mut_count = {}
    mut_id = {}
    for id_ in id_list:
        try:
            # select column for this id
            column = vcf_df[id_]
            # go through all values for this id
            for row,value in enumerate(column):
                # if value does not start with point there is a mutation
                if not value.startswith("."):
                    # if there is muation we save the row and count
                    if row in mut_count.keys():
                        mut_count[row] += 1
                        mut_id[row].append(id_)
                    else:
                        mut_count[row] = 1
                        mut_id[row] = [id_]
                    
        except KeyError:
            print(f"{id_} not in vcf")
    
    mut_count = {vcf_df.POS[k]:v for k,v in mut_count.items()}
    mut_id = {vcf_df.POS[k]:v for k,v in mut_id.items()}
            
    return mut_count, mut_id

def identify_gene(pos_dic,gff_df):
    for pos in pos_dic.keys():
        
    
    return



mut_count_1, mut_id_1 = get_mut(cluster_1_id.new_names, vcf)
mut_count_2, mut_id_2 = get_mut(cluster_2_id.new_names, vcf)

# check if all mutations appear atlest once -> only from parent not
all_mut = []
for mut in mut_count_1.keys():
    if mut not in all_mut:
        all_mut.append(mut)
        
for mut in mut_count_2.keys():
    if mut not in all_mut:
        all_mut.append(mut)   
        
for mut in vcf.POS:
    if mut not in all_mut:
        print(f"{mut} not in all mut")


# import gff
gff = read_gff("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff")

# look for the positions in gff file and assign a gene to the vcf position
genes_1 = identify_gene(mut_count_1, gff)




