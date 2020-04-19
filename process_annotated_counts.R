library("data.table", quietly=TRUE)
library('stringr', quietly=TRUE)
library("dplyr", quietly=TRUE)
library("methylKit", quietly=TRUE)

p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
get_sample_mean <- function(count_table, epi_id1, epi_id2){
  print("---> Calculating row means")
  count_table = data.frame(rowMeans(count_table[,grep(epi_id1, colnames(count_table))]),
                           rowMeans(count_table[,grep(epi_id2, colnames(count_table))]))
  colnames(count_table) = c(epi_id1, epi_id2)
  return(count_table)
}
get_cor_sample_pair <- function(count_table, epi_id1, epi_id2){
  print("---> Calculating pairwise correlation")
  count_table.cor <- count_table %>%
    mutate(gene = str_split_fixed(rownames(count_table), "\\:", 2)[1:dim(count_table)[[1]]]) %>%
    group_by(gene) %>%
    mutate(freq = length(gene)) %>%
    filter(freq >= 3)  %>%
    mutate(r_val = cor(!!as.name(epi_id1), !!as.name(epi_id2))) %>%
    filter(!is.na(r_val)) %>%
    dplyr::select(gene, freq, r_val) %>%
    unique() %>%
    mutate(p_val = p_value_calculator(r_val, freq),
           p_adj = p.adjust(p_val, method="fdr", n=freq))
  return(count_table.cor)
}
summarize_cor <- function(cor_table, padj, rval){
  anticorrelating_genes = c(cor_table[cor_table$`p_val` <= padj & cor_table$`r_val` <= -rval,]$`gene`)
  noncorrelating_genes = c(cor_table[cor_table$`p_val` <= padj & cor_table$`r_val` < rval & cor_table$`r_val` > -rval,]$`gene`)
  correlating_genes = c(cor_table[cor_table$`p_val` <= padj & cor_table$`r_val` >= rval,]$`gene`)
  return(list(anticorrelating_genes, noncorrelating_genes, correlating_genes))
}
prepare_methylation_count <- function(methylation.table){
  methylation.table = data.frame(methylation.table[,c(1,4,5,7,9,14,15)])
  methylation.table$V14 = as.numeric(as.character(methylation.table$V14))
  methylation.table$V15 = as.numeric(as.character(methylation.table$V15))
  methylation.table$V14[is.na(methylation.table$V14)] = 0
  methylation.table$V15[methylation.table$V15 < 0] = 0
  methylation.table$V15[is.na(methylation.table$V15)] = 0
  methylation.table.aggregate <- aggregate(list(readcov_1 = methylation.table$V14, readcov_2 = methylation.table$V15),
                                           FUN=mean,
                                           by = list(chr = methylation.table$V1, start = methylation.table$V4,
                                                     end = methylation.table$V5, strand = methylation.table$V7, gene_id = methylation.table$V9))
  new_id = as.data.table(str_split_fixed(methylation.table.aggregate$gene_id, "\\;|\\ ", 8))[,c(2,8)]
  new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
  new_id = data.table(paste(new_id$V2,new_id$V8, sep=':'))
  
  count = data.table(new_id$V1, methylation.table.aggregate$readcov_1, methylation.table.aggregate$readcov_2)
  count = data.frame(count[order(count$V1),])
  rownames(count) = count$V1
  count = count[,2:3]
  colnames(count) = c(epi_id1, epi_id2)
  return(count)
}
prepare_histone_count <- function(histone.table){
  histone.table = data.frame(histone.table[,c(1,4,5,7,9,13,14)])
  histone.table$V13 = as.numeric(as.character(histone.table$V13))
  histone.table$V14 = as.numeric(as.character(histone.table$V14))
  histone.table$V13[is.na(histone.table$V13)] = 0
  histone.table$V14[histone.table$V14 < 0] = 0
  histone.table.aggregate <- aggregate(list(readcov_1 = histone.table$V13, readcov_2 = histone.table$V14),
                                       FUN=mean,
                                       by = list(chr = histone.table$V1, start = histone.table$V4,
                                                 end = histone.table$V5, strand = histone.table$V7, gene_id = histone.table$V9))
  new_id = as.data.table(str_split_fixed(histone.table.aggregate$gene_id, "\\;|\\ ", 8))[,c(2,8)]
  new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
  new_id = data.table(paste(new_id$V2,new_id$V8, sep=':'))
  
  count = data.table(new_id$V1, histone.table.aggregate$readcov_1, histone.table.aggregate$readcov_2)
  count = data.frame(count[order(count$V1),])
  rownames(count) = count$V1
  count = count[,2:3]
  # colnames(count) = c(epi_id1, epi_id2)
  return(count)
}

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine")
gene_id = data.frame(fread("gene_id", header=FALSE, col.names=c("gene_id")))
epi_id1 = "E003"
epi_id2 = "E004"
epi_id3 = "E005"

#================PREPARE COUNTFILES===============
prepare_countfile <- function(feature, epi_id1, epi_id2){
  file_path = normalizePath(paste(getwd(), feature, "annotatedcounts", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  
  if (feature == "expression") {
    expression.countfile = list.files(file_path, pattern=pair, full.names=TRUE)
    print("---> Reading count table")
    expression.count = data.frame(fread(expression.countfile[[1]], header=TRUE,stringsAsFactors=FALSE))
    rownames(expression.count) = gene_id$gene_id
    expression.count = get_sample_mean(expression.count, epi_id1, epi_id2)
    print("---> Writing count table into file")
    fwrite(expression.count,
           file = paste(file_path, paste(epi_id1, epi_id2, "counts.txt", sep='_'), sep = '/'),
           row.names = TRUE, col.names = TRUE, sep='\t', quote=FALSE)
    return(expression.count)
  }
  else if (feature == "methylation"){
    methylation.countfile = list.files(file_path, pattern=pair, full.names=TRUE)
    print("---> Reading count table")
    methylation.table = fread(methylation.countfile[[1]], header=FALSE, stringsAsFactors = FALSE)
    print("---> Compressing count table")
    methylation.count = prepare_methylation_count(methylation.table)
    colnames(methylation.count) = c(epi_id1, epi_id2)
    print("---> Writing count table into file")
    fwrite(methylation.count,
           file = paste(file_path, paste(epi_id1, epi_id2, "counts.txt", sep='_'), sep = '/'),
           row.names = TRUE, col.names = TRUE, sep='\t', quote=FALSE)
    return(methylation.count)
  }
  else {
    histone.countfile = list.files(file_path, pattern=pair, full.names=TRUE)
    print("---> Reading count table")
    histon.table = fread(histone.countfile[[1]], header=FALSE, stringsAsFactors = FALSE)
    print("---> Compressing count table")
    histone.count = prepare_histone_count(histon.table)
    colnames(histone.count) = c(epi_id1, epi_id2)
    print("---> Writing count table into file")
    fwrite(histone.count,
           file = paste(file_path, paste(epi_id1, epi_id2, "counts.txt", sep='_'), sep = '/'),
           row.names = TRUE, col.names = TRUE, sep='\t', quote=FALSE)
    return(histone.count)
  }
}
expression.count1 = prepare_countfile("expression", epi_id1, epi_id2)
expression.count2 = prepare_countfile("expression", epi_id1, epi_id3)
expression.count3 = prepare_countfile("expression", epi_id2, epi_id3)

methylation.count1 = prepare_countfile("methylation", epi_id1, epi_id2)
methylation.count2 = prepare_countfile("methylation", epi_id1, epi_id3)
methylation.count3 = prepare_countfile("methylation", epi_id2, epi_id3)

H3K36me3.count1 = prepare_countfile("H3K36me3", epi_id1, epi_id2)
H3K36me3.count2 = prepare_countfile("H3K36me3", epi_id1, epi_id3)
H3K36me3.count3 = prepare_countfile("H3K36me3", epi_id2, epi_id3)

H3K27ac.count1 = prepare_countfile("H3K27ac", epi_id1, epi_id2)
H3K27ac.count2 = prepare_countfile("H3K27ac", epi_id1, epi_id3)
H3K27ac.count3 = prepare_countfile("H3K27ac", epi_id2, epi_id3)


#===========================================
#1: ca 2 >= 0.75, 2: 1 > 0.75, 3: ca 2 duoi 0.75
methylKit:::readMethylBaseDB("/Users/dhthutrang/metDB/methylDB\ 2020-04-06\ M3E/methylDiff_E003_E005.txt.bgz",
                             sample.ids = c("E003", "E004"), assembly = 'hg19', context="CpG", resolution="base",
                             dbtype="tabix", treatment = c(0,1), destranded=FALSE)


