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
from collections import defaultdict
import gzip
import re
import itertools
import igraph as ig
import matplotlib.pyplot as plt

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



GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')
                      
def read_gff(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)

def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

    
                      


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


def mutation_graph(mut_id_1, mut_id_2, df):
    
    list_clones = list(df.new_names)
    
    n_vertices = len(list_clones)
    edges = []
    edge_color = []
    
    total_dic = {}
    # combine both dictionaries
    for key_1 in mut_id_1.keys():
        if key_1 in mut_id_2.keys():
            total_dic[key_1] = mut_id_1[key_1] + mut_id_2[key_1]
        else:
            total_dic[key_1] = mut_id_1[key_1]
    for key_2 in mut_id_2.keys():
        if key_2 in total_dic.keys():
            continue
        else:
            total_dic[key_2] = mut_id_2[key_2]
    
    # generate edges list
    for clones in list(total_dic.values()):
        # get index in clones list
        c = [list_clones.index(clone) for clone in clones]
        # generate all egdes
        combinations = list(itertools.combinations(c,2))
        if len(combinations) == 0:
            combinations = [(c[0],c[0])]
        # add to edges list
        edges.extend(combinations)
        
    for pair in edges:
        vert_1 = pair[0]
        vert_2 = pair[1]
        
        cond_1 = df["Condition"][vert_1]
        cond_2 = df["Condition"][vert_2]
        
        if cond_1 == cond_2:
            if cond_1 == "supernatant":
                edge_color.append("orange")
            else:
                edge_color.append("blue")
        else:
            edge_color.append("red")
    
    # create igraph
    g = ig.Graph(n_vertices, edges)
    
    # color edges
    g.es["color"] = edge_color
    
    # Set attributes for the graph, nodes, and edges
    g["title"] = "Mutation Network"
    g.vs["name"] = list_clones
    g.vs["cluster"] = df["phenotype_clusters"]
    # g.vs["x"] = df["PC1"]
    # g.vs["y"] = df["PC2"]
    # todo: add further attribures as "unique", "pathway"
    g.vs["condition"] = df["Condition"]
    
    
    return g

def visualize_graph(g,coord):
    fig, ax = plt.subplots(figsize=(5,5))
    ig.plot(
    g,
    layout="circle",
    target=ax,
    # layout=coord,
    vertex_size=0.05,
    vertex_color=["steelblue" if cluster == 1 else "salmon" for cluster in g.vs["cluster"]],
    vertex_frame_width=2.0,
    vertex_frame_color="white",
    vertex_label=g.vs["name"],
    vertex_label_size=7.0,
    vertex_label_dist = 10,
    vertex_label_color = "blue",
    edge_width = 1,
    # edge_width=[2 if married else 1 for married in g.es["married"]],
    
    )

    plt.show()
    
    return fig
                

def identify_function(genes_dic,gff_df):
    
    functions_dic = {}
    
    for pos,value in genes_dic.items():
        count = value[0]
        ids = value[1]

        function = gff_df["product"][pos+1]
        
        if function not in functions_dic.keys():
            functions_dic[function] = [count, [pos], ids]
        else:
            functions_dic[function][0] += count
            functions_dic[function][1].append(pos)
            for id_ in ids.keys():
                if id_ not in functions_dic[function][2].keys():
                    functions_dic[function][2][id_] = 1
                else:
                    functions_dic[function][2][id_] += ids[id_]        
    
    return functions_dic


# import vcf
vcf = read_vcf("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/2020-09-25_mod_SA6850_filteredFromAncestor.vcf")                     
# import gff
gff = read_gff("/Users/bp/Uni/Computational/HS22/BIO253/Data/SA6850_set_forIGV/SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff")
# import clustering results
df = pd.read_csv('/Users/bp/Uni/Computational/HS22/BIO253/Bio253_research_cycle_genomics/Out/SA6850/cluster_results_SA6850.csv')
# adjust clones id
df = change_name(df)


# choose list of id for each cluster
cluster_1_id = df[df["phenotype_clusters"] == 1]
cluster_2_id = df[df["phenotype_clusters"] == 2]

# get the mutations
mut_count_1, mut_id_1 = get_mut(cluster_1_id.new_names, vcf)
mut_count_2, mut_id_2 = get_mut(cluster_2_id.new_names, vcf)
# todo: get rid of synonymous mutations


# look for the positions in gff file and assign a gene to the vcf position
genes_1 = identify_gene(mut_count_1, mut_id_1,gff)
genes_2 = identify_gene(mut_count_2, mut_id_2,gff)


# visualize mutations shared
genes_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys()}
genes_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys()}
g = mutation_graph(genes_1_id, genes_2_id, df)
coord = np.array(list(zip(df["PC1"],df["PC2"])))
fig = visualize_graph(g,coord)
fig.savefig('../Out/SA6850/phenotype_SA6850_all.eps', format='eps')

# get symmetric diffences in genes
unique_1 = {row:(genes_1[row]) for row in genes_1.keys() if row not in genes_2.keys()}
unique_2 = {row:(genes_2[row]) for row in genes_2.keys() if row not in genes_1.keys()}

# product level
unique_function_1 = identify_function(unique_1, gff)
unique_function_2 = identify_function(unique_2, gff)


# visualize only unique mutations
unique_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys() if row not in genes_2.keys()}
unique_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys() if row not in genes_1.keys()}
g_unique = mutation_graph(unique_1_id, unique_2_id, df)
fig = visualize_graph(g_unique,coord)
fig.savefig('../Out/SA6850/phenotype_SA6850_unique.eps', format='eps')



# visualize functional level
unique_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys() if row not in genes_2.keys()}
unique_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys() if row not in genes_1.keys()}
g_unique = mutation_graph(unique_1_id, unique_2_id, df)
fig = visualize_graph(g_unique,coord)
fig.savefig('../Out/SA6850/phenotype_SA6850_function_unique.eps', format='eps')


#todo: get unique combination of mutations


# todo: add analytical measurements for graph


# get common mutations
common = {row:[(genes_1[row][0] , genes_2[row][0]), list(genes_1[row][1].keys()) + list(genes_2[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys()}


# sort dictionaries we want to print
genes_1 = sorted(genes_1.items(), key=lambda x:x[0])
genes_2 = sorted(genes_2.items(), key=lambda x:x[0])
unique_1 = sorted(unique_1.items(), key=lambda x:x[0])
unique_2 = sorted(unique_2.items(), key=lambda x:x[0])
common = sorted(common.items(), key=lambda x:x[0])


# create ranges
ranges = {}
for index, row in gff.iterrows():
    if index == 0:
        continue
    # only exract the "gene" row
    if f"{row.start}-{row.end}" not in ranges.keys() and row.feature == "gene":
        ranges[index] = f"{row.start}-{row.end}"

# write it to file
with open("../Out/SA6850/mutated_genes_SA6850_phenotype.txt", 'w') as f:
    f.write(f"\n#Cluster 1 ({len(genes_1)}):\n")
    for key, value in genes_1: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    f.write(f"\n\n#Cluster 2 ({len(genes_2)}):\n")
    for key, value in genes_2: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    f.write(f"\n\n#Unique to Cluster 1 ({len(unique_1)}):\n")
    for key,value in unique_1:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    f.write(f"\n\n#Unique to Cluster 2 ({len(unique_2)}):\n")
    for key,value in unique_2:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    f.write(f"\n\n#Common ({len(common)}):\n")
    for key,value in common:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
        
# store gff file
gff.to_csv('../Out/gff_SA6850.csv')


def extract_condition(dic):
    conditions = []
    for id_ in dic.keys():
        cond = id_.split("_")[1]
        if cond == "S":
            conditions.append("S")
        else:
            conditions.append("T")
    return conditions
        


# store mutations list in csv
with open("../Out/SA6850/list_mutated_genes_SA6850_phenotype.csv", 'w') as f:
    f.write("Cluster,Mutated_gene,Condition\n")
    for key, value in genes_1: 
        f.write('1,%s,T\n' % (key))
    for key, value in genes_2: 
        cond = extract_condition(value[1])
        if "S" in cond and "T" in cond:
            f.write('2,%s,T/S\n' % (key))
        elif "S" in cond:
            f.write('2,%s,S\n' % (key))
        else:
            f.write('2,%s,T\n' % (key))


    