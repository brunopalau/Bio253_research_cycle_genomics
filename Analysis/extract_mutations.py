#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:09:54 2022

@author: bp
"""

import pandas as pd
import numpy as np
import json

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

def identify_gene(pos_dic,id_dic,gff_df):
    # create dict of ranges:`rows assigned to the range`
    ranges = {}
    for index, row in gff_df.iterrows():
        if index == 0:
            continue
        # only exract the "gene" row
        if f"{row.start}-{row.end}" not in ranges.keys() and row.feature == "gene":
            ranges[f"{row.start}-{row.end}"] = index
    
    # go through mutation and get rows in gff to which it belongs to
    genes_dic = {}
    for pos, count in pos_dic.items():
        # go through all possible ranges
        for ran,row in ranges.items():
            start,stop = ran.split("-")
            # if position of mutation is part of the gene than add {gene:(count,list_clones)}
            if pos <= int(stop) and pos >= int(start):
                # store only first row of a unique range
                if row not in genes_dic.keys():
                    genes_dic[row] = [count,{key:1 for key in id_dic[pos]}]
                else:
                    genes_dic[row][0] += count
                    for id_ in id_dic[pos]:
                        if id_ in genes_dic[row][1].keys():
                            genes_dic[row][1][id_] += 1
                        else:
                            genes_dic[row][1][id_] = 1
    return genes_dic

# import vcf
vcf = read_vcf("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/2020-09-25_mod_SA6850_filteredFromAncestor.vcf")                     
# import gff
gff = read_gff("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff")
# import clustering results
df = pd.read_csv('/Users/bp/Uni/Computational/HS22/BIO253/Data/cluster_results.csv')
# adjust clones id
df = change_name(df)


# choose list of id for each cluster
cluster_1_id = df[df["phenotype_clusters_2"] == 1]
cluster_2_id = df[df["phenotype_clusters_2"] == 2]

# get the mutations
mut_count_1, mut_id_1 = get_mut(cluster_1_id.new_names, vcf)
mut_count_2, mut_id_2 = get_mut(cluster_2_id.new_names, vcf)
# todo: get rid of synonymous mutations


# look for the positions in gff file and assign a gene to the vcf position
genes_1 = identify_gene(mut_count_1, mut_id_1,gff)
genes_2 = identify_gene(mut_count_2, mut_id_2,gff)

# get symmetric diffences in genes
unique_1 = {row:(genes_1[row]) for row in genes_1.keys() if row not in genes_2.keys()}
unique_2 = {row:(genes_2[row]) for row in genes_2.keys() if row not in genes_1.keys()}

# get unique combination of mutations



# get common mutations
common = {row:[genes_1[row][0] + genes_1[row][0], list(genes_1[row][1].keys()) + list(genes_1[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys()}


# sort dictionaries we want to print
genes_1 = sorted(genes_1.items(), key=lambda x:x[0])
genes_2 = sorted(genes_2.items(), key=lambda x:x[0])
unique_1 = sorted(unique_1.items(), key=lambda x:x[0])
unique_2 = sorted(unique_2.items(), key=lambda x:x[0])
common = sorted(common.items(), key=lambda x:x[0])


# write it to file
with open("../Out/mutated_genes_SA6850_phenotypes.txt", 'w') as f:
    f.write(f"\n#Cluster 1 ({len(genes_1)}):\n")
    for key, value in genes_1: 
        f.write('%s:%s\n' % (key, value))
    f.write(f"\n\n#Cluster 2 ({len(genes_2)}):\n")
    for key, value in genes_2: 
        f.write('%s:%s\n' % (key, value))
    f.write(f"\n\n#Unique to Cluster 1 ({len(unique_1)}):\n")
    for key,value in unique_1:
        f.write('%s:%s\n' % (key, value))
    f.write(f"\n\n#Unique to Cluster 2 ({len(unique_2)}):\n")
    for key,value in unique_2:
        f.write('%s:%s\n' % (key, value))
    f.write(f"\n\n#Common ({len(common)}):\n")
    for key,value in common:
        f.write('%s:%s\n' % (key, value))
        
        
        
        
        
        