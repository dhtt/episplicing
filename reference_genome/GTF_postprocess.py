#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:10:44 2019

@author: dohoangthutrang
"""
import pandas as pd
import numpy as np
import HTSeq
import itertools
import sys
import json, pickle, collections
#
gtf_file = sys.argv[1]
out_file = sys.argv[2]
gene_id_file = sys.argv[3]

# Get table containing trancript and gene_id
with open(gene_id_file) as f:
    coding_gene = pd.read_table(f, sep= '\t', index_col=None, header=0, lineterminator='\n')
    coding_gene.columns = ['transcript_id', 'gene_id', 'description']


print("=====> STRIPPING SPACES FROM GENE ID")
for i in range(len(coding_gene['gene_id'])):
    if (' ' in coding_gene['gene_id'][i]):
        print(i, coding_gene['gene_id'][i])
        coding_gene['gene_id'][i] = coding_gene['gene_id'][i].replace(' ','')
        print(i, coding_gene['gene_id'][i])
        
# Convert transcipt to gene_id
print("=====> CONVERTING TRANSCRIPT TO GENE ID")
gene_id = {} #store number of time a gene appears (key:gene, val: counts)
transcipttogenes_id = {} #for faster lookup of existed rtanscript (key: original transcript to be converted, val: associated gene_id) 
all_feat = []
for f in HTSeq.GFF_Reader( gtf_file ):
    if f.type == "aggregate_gene": 
        genes = [] #store converted transcript as a list of gene ids to be joined later into string 
        transcript_id = f.attr['gene_id']
        for transcript in (f.attr['gene_id'].split("+")):
            if coding_gene[coding_gene['transcript_id'] == transcript]['gene_id'].iloc[0] not in genes: #avoid cases like SAMD1+SAMD1+....
                genes.append(coding_gene[coding_gene['transcript_id'] == transcript]['gene_id'].iloc[0])
        genes = "+".join(genes)
    
        if genes not in gene_id.keys(): #if gene appears once
            f.attr['gene_id'] = genes
            gene_id[genes] = 1
        else: #if gene appears more than once
            f.attr['gene_id'] = "__00".join([genes, str(gene_id[genes])])
            gene_id[genes] += 1
        transcipttogenes_id[transcript_id] = f.attr['gene_id']
    else:
        f.attr['gene_id'] = transcipttogenes_id[transcript_id]
    f.name = f.attr['gene_id']
    
    #Store f as GenomicFeature with new gene_id
    feat = HTSeq.GenomicFeature( f.name, f.type, 
                                HTSeq.GenomicInterval( f.iv.chrom, f.iv.start, f.iv.end, f.iv.strand ) )
    feat.attr = {}
    feat.attr['gene_id'] = f.attr['gene_id']
    if f.type == "exonic_part":
        feat.attr['transcripts'] = transcript_id
        feat.attr['exonic_part_number'] = f.attr['exonic_part_number']
    all_feat.append(feat)
    
###Step 5: Sort the aggregates, then write everything out
print("=====> WRITING RESULT TO FILE")
fout = open( out_file, "w" ) 
for feat in all_feat:
   fout.write( feat.get_gff_line() )
fout.close()    