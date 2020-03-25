library("tidyr")
library("genomation")
library("data.table")
library('NMF')
library('ggplot2')
library('gplots')
library('reshape')
library('Hmisc')
library('stringr')
library("RColorBrewer")
library("dplyr")
library("DEXSeq")
suppressPackageStartupMessages( library( "DEXSeq" ) )
setwd("~/Documents/BIOINFO/Episplicing/files/Result/mRNA/rep3")

normalize_01<-function(m){ 
  (m - min(m))/(max(m)-min(m))
}
normalize_11<-function(m){
  2*(m - min(m))/(max(m)-min(m)) - 1
}

new_label = c("H1 Cell Line (ESC H1)", "iPSc 19", "iPSc 6","Mesenchymal Stem Cell", "Mesendoderm", "Trophoblast", "Neural Progenitor Cell", "IMR90","Fetal Brain", "Fetal Brain Germinal Matrix", "Brain Hippocampus Middle", "Breast Luminal Epithelium",  "Breast Myoepithelial cells", "Adult Liver")
gene_list_anticor = c("ENY2", "ADAM9", "BNIP3L", "LIPG", "SLC16A2", "MRPS18C", "MRPS30", "SLC25A43", "FOXO4", "^CA2")
gene_list_cor = c("MRPL22", "PQBP1", "SLC1A3", "NSG1", "PRR7", "TMED9", "RHPN1", "APCDD1", "UGDH", "FOXO4", "^CA2")
cell_list_grep = paste(cell_types, collapse = "|")
gene_list_anticor_grep = paste(gene_list_anticor, collapse = "|")
gene_list_cor_grep = paste(gene_list_cor, collapse = "|")


inDir = normalizePath(paste(getwd(),"count", sep="/"))
countFiles = list.files(inDir, pattern="*count.txt", full.names=TRUE)
basename(countFiles)
file_name = as.data.table(str_split_fixed(basename(countFiles), "\\.", 5))
row_names = c(str_split_fixed(basename(countFiles), "\\.bed", 5)[,1])
condition = file_name$V2
cell_types = c("H1hesc","Huvec")

sampleTable = data.frame(
  row.names = row_names,
  condition = condition)

dxd_stranded = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile= normalizePath(paste(getwd(),"reference.flattened.ens.gtf", sep="/"))
)
head(counts(dxd_stranded))
dxd_stranded.e1 <- estimateSizeFactors(dxd_stranded)
dxd_stranded.e2 = estimateDispersions(dxd_stranded.e1)
dxd_stranded.deu = testForDEU(dxd_stranded.e2)
dxd_stranded.fc = estimateExonFoldChanges( dxd_stranded.deu, fitExpToVar="condition")
dxd_stranded.res = DEXSeqResults(dxd_stranded.fc)

table ( dxr1$padj < 0.1 )
table (tapply( dxr1$padj < 0.1, dxr1$groupID, any ))