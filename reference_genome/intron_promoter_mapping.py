#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 09:37:53 2019

@author: dohoangthutrang
"""
import sys, collections, itertools, os.path
import HTSeq

#gtf_file = 'UCSC_coding_main_ens_flattened_temp.gtf'

gtf_file = sys.argv[1]
promoter_file = sys.argv[2]
intron_file = sys.argv[3]

intron_part = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
promoter_part = HTSeq.GenomicArrayOfSets( "auto", stranded=True )
for f in HTSeq.GFF_Reader(gtf_file):
    f.attr['gene_id'] = f.attr['gene_id'].replace( ":", "_" )
    if f.type != 'aggregate_gene':
        intron_part[f.iv] += ( f.attr['gene_id'])
    else:
        promoter_part[f.iv] += ( f.attr['gene_id'])
     
#Annotate promoters
promoters = []
for iv, s in promoter_part.steps():
    if len(s) == 0:
      continue
    if iv.strand == "+":
        new_iv = HTSeq.GenomicInterval( iv.chrom, iv.start-2000, iv.start, iv.strand )
    else:
        new_iv = HTSeq.GenomicInterval( iv.chrom, iv.end, iv.end+2000, iv.strand )
    
    for g_id in s:
        gene_id = g_id
    promoter = HTSeq.GenomicFeature( gene_id, "promoter", new_iv)
    promoter.attr = {'gene_id': gene_id}
    promoters.append(promoter)
promoters.sort( key = lambda promoter: ( promoter.iv.chrom, promoter.iv.start ) )
fout = open(promoter_file, "w" ) 
for promoter in promoters:
   fout.write(promoter.get_gff_line())
fout.close()

#Annotate introns
introns = []
for iv, s in intron_part.steps():
    if len(s) == 0:
      continue
    iv = HTSeq.GenomicInterval( iv.chrom, iv.start, iv.end, iv.strand )
    for g_id in s:
        gene_id = g_id
    intron = HTSeq.GenomicFeature( gene_id, "promoter", iv)
    intron.attr = {'gene_id': gene_id}
    introns.append(intron)

intron_list = []
i = 0 #Counter for intronic_parts that belong to same gene cluster
for m in range(len(introns)-1):
    if introns[m].name == introns[m+1].name:
        new_iv = HTSeq.GenomicInterval( introns[m].iv.chrom, introns[m].iv.end, introns[m+1].iv.start, introns[m+1].iv.strand )
        intron = HTSeq.GenomicFeature( introns[m].name, "intron", new_iv)
        intron.attr = {'gene_id': introns[m].name, 'intronic_part_number': "%03d" % ( i+1 )}
        intron_list.append(intron)
        i+=1
    else: 
        i = 0
intron_list.sort( key = lambda intron: ( intron.iv.chrom, intron.iv.start ) )
fout = open(intron_file, "w" ) 
for intron in intron_list:
   fout.write(intron.get_gff_line())
fout.close()