#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 20:30:28 2022

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
        
        new_string = "JE2_"
        
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


def mutation_graph(mut_id_1, mut_id_2, mut_id_3, mut_id_4, df):
    
    list_clones = list(df.new_names)
    
    n_vertices = len(list_clones)
    edges = []
    edge_color = []
    
    total_dic = {}
    
    # list with all pos
    positions = list(mut_id_1) + list(mut_id_2) + list(mut_id_3) + list(mut_id_4) 
    
    # create list of shared mutations
    for pos in positions:
        if pos not in total_dic.keys():
            total_dic[pos] = []
        
        try:
            total_dic[pos] += mut_id_1[pos]
        except: pass
        try:
            total_dic[pos] += mut_id_2[pos]
        except: pass
        try:
            total_dic[pos] += mut_id_3[pos]
        except: pass
        try:
            total_dic[pos] += mut_id_4[pos]
        except: pass

    
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
        
        
    # color dictionary
    # color_dict = {(1,2): "blue", (2,1): "blue", (1,3): "orange", (3,1): "orange", (1)}
        
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
    vertex_color=["steelblue" if cluster == 1 else "red" if cluster == 2 else "salmon" if cluster ==3 else "black" for cluster in g.vs["cluster"]],
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
vcf = read_vcf("/Users/bp/Uni/Computational/HS22/BIO253/Data/JE2_set_forIGV/2020-09-25_mod_SAJE2_filteredFromAncestor_noOutlier.1.vcf")                     
# import gff
gff = read_gff("/Users/bp/Uni/Computational/HS22/BIO253/Data/JE2_set_forIGV/SE_SE2_GCA_002085525.1_ASM208552v1_genomic.gff")
# import clustering results
df = pd.read_csv('/Users/bp/Uni/Computational/HS22/BIO253/Bio253_research_cycle_genomics/Out/JE2/cluster_results_JE2.csv')
# adjust clones id
df = change_name(df)


# choose list of id for each cluster
cluster_1_id = df[df["phenotype_clusters"] == 1]
cluster_2_id = df[df["phenotype_clusters"] == 2]
cluster_3_id = df[df["phenotype_clusters"] == 3]
cluster_4_id = df[df["phenotype_clusters"] == 4]


# get the mutations
mut_count_1, mut_id_1 = get_mut(cluster_1_id.new_names, vcf)
mut_count_2, mut_id_2 = get_mut(cluster_2_id.new_names, vcf)
mut_count_3, mut_id_3 = get_mut(cluster_3_id.new_names, vcf)
mut_count_4, mut_id_4 = get_mut(cluster_4_id.new_names, vcf)
# todo: get rid of synonymous mutations


# look for the positions in gff file and assign a gene to the vcf position
genes_1 = identify_gene(mut_count_1, mut_id_1,gff)
genes_2 = identify_gene(mut_count_2, mut_id_2,gff)
genes_3 = identify_gene(mut_count_3, mut_id_3,gff)
genes_4 = identify_gene(mut_count_4, mut_id_4,gff)


# visualize mutations shared
genes_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys()}
genes_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys()}
genes_3_id = {row:(list(genes_3[row][1].keys())) for row in genes_3.keys()}
genes_4_id = {row:(list(genes_4[row][1].keys())) for row in genes_4.keys()}

# create and plot graph
g = mutation_graph(genes_1_id, genes_2_id, genes_3_id, genes_4_id,df)
coord = np.array(list(zip(df["PC1"],df["PC2"])))
fig = visualize_graph(g,coord)
fig.savefig('../Out/JE2/phenotype_JE2_all.eps', format='eps')

# get symmetric diffences in genes
unique_1 = {row:(genes_1[row]) for row in genes_1.keys() if row not in genes_2.keys() and row not in genes_3.keys() and row not in genes_4.keys()}
unique_2 = {row:(genes_2[row]) for row in genes_2.keys() if row not in genes_1.keys() and row not in genes_3.keys() and row not in genes_4.keys()}
unique_3 = {row:(genes_3[row]) for row in genes_3.keys() if row not in genes_1.keys() and row not in genes_2.keys() and row not in genes_4.keys()}
unique_4 = {row:(genes_4[row]) for row in genes_4.keys() if row not in genes_1.keys() and row not in genes_2.keys() and row not in genes_3.keys()}

# product level
unique_function_1 = identify_function(unique_1, gff)
unique_function_2 = identify_function(unique_2, gff)
unique_function_3 = identify_function(unique_3, gff)
unique_function_4 = identify_function(unique_4, gff)



# visualize only unique mutations
unique_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys() if row not in genes_2.keys() and row not in genes_3.keys() and row not in genes_4.keys()}
unique_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys() if row not in genes_1.keys() and row not in genes_3.keys() and row not in genes_4.keys()}
unique_3_id = {row:(list(genes_3[row][1].keys()))for row in genes_3.keys() if row not in genes_1.keys() and row not in genes_2.keys() and row not in genes_4.keys()}
unique_4_id = {row:(list(genes_4[row][1].keys()))for row in genes_4.keys() if row not in genes_1.keys() and row not in genes_2.keys() and row not in genes_3.keys()}

g_unique = mutation_graph(unique_1_id, unique_2_id, unique_3_id, unique_4_id, df)
fig = visualize_graph(g_unique,coord)
fig.savefig('../Out/JE2/phenotype_JE2_unique.eps', format='eps')



# visualize functional level
# unique_1_id = {row:(list(genes_1[row][1].keys())) for row in genes_1.keys() if row not in genes_2.keys()}
# unique_2_id = {row:(list(genes_2[row][1].keys())) for row in genes_2.keys() if row not in genes_1.keys()}
# g_unique = mutation_graph(unique_1_id, unique_2_id, unique_3_id, unique_4_id, df)
# visualize_graph(g_unique,coord)

#todo: get unique combination of mutations


# todo: add analytical measurements for graph


# get common mutations
common = {row:[(genes_1[row][0] , genes_2[row][0]), list(genes_1[row][1].keys()) + list(genes_2[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys() and row in genes_3.keys() and row in genes_4.keys()}
# common = {row:[(genes_1[row][0] , genes_2[row][0]), list(genes_1[row][1].keys()) + list(genes_2[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys() and row in genes_3.keys() and row in genes_4.keys()}
# common = {row:[(genes_1[row][0] , genes_2[row][0]), list(genes_1[row][1].keys()) + list(genes_2[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys() and row in genes_3.keys() and row in genes_4.keys()}
# common = {row:[(genes_1[row][0] , genes_2[row][0]), list(genes_1[row][1].keys()) + list(genes_2[row][1].keys())] for row in genes_1.keys() if row in genes_2.keys() and row in genes_3.keys() and row in genes_4.keys()}


# sort dictionaries we want to print
genes_1 = sorted(genes_1.items(), key=lambda x:x[0])
genes_2 = sorted(genes_2.items(), key=lambda x:x[0])
genes_3 = sorted(genes_3.items(), key=lambda x:x[0])
genes_4 = sorted(genes_4.items(), key=lambda x:x[0])
unique_1 = sorted(unique_1.items(), key=lambda x:x[0])
unique_2 = sorted(unique_2.items(), key=lambda x:x[0])
unique_3 = sorted(unique_3.items(), key=lambda x:x[0])
unique_4 = sorted(unique_4.items(), key=lambda x:x[0])
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
with open("../Out/JE2/mutated_genes_JE2_phenotype.txt", 'w') as f:
    # cluster 1
    f.write(f"\n#Cluster 1 ({len(genes_1)}):\n")
    for key, value in genes_1: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    # cluster 2
    f.write(f"\n\n#Cluster 2 ({len(genes_2)}):\n")
    for key, value in genes_2: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    # cluster 3
    f.write(f"\n\n#Cluster 3 ({len(genes_3)}):\n")
    for key, value in genes_3: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    # cluster 4
    f.write(f"\n\n#Cluster 4 ({len(genes_4)}):\n")
    for key, value in genes_4: 
        f.write('%s,%s,%s:%s\n' % (key,gff.gene[key],ranges[key], value))
    # unique cluster 1
    f.write(f"\n\n#Unique to Cluster 1 ({len(unique_1)}):\n")
    for key,value in unique_1:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    # unique cluster 2
    f.write(f"\n\n#Unique to Cluster 2 ({len(unique_2)}):\n")
    for key,value in unique_2:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    # unique cluster 3
    f.write(f"\n\n#Unique to Cluster 3 ({len(unique_3)}):\n")
    for key,value in unique_3:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    # unique cluster 4
    f.write(f"\n\n#Unique to Cluster 4 ({len(unique_4)}):\n")
    for key,value in unique_4:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))
    # common
    f.write(f"\n\n#Common ({len(common)}):\n")
    for key,value in common:
        f.write('%s,%s,%s: %s\n' % (key,gff.gene[key],ranges[key], value))

# store gff file
gff.to_csv('../Out/JE2/gff_JE2.csv')
