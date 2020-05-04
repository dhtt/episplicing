library("data.table", quietly=TRUE)
library('stringr', quietly=TRUE)
library("plyr", quietly=TRUE)
library("dplyr", quietly=TRUE)
library(GenomicRanges)

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
  histone.table = data.frame(histone.table[,c(1,4,5,7,9,15,16)])
  histone.table$V15 = as.numeric(as.character(histone.table$V15))
  histone.table$V16 = as.numeric(as.character(histone.table$V16))
  histone.table$V15[is.na(histone.table$V15)] = 0
  histone.table$V16[histone.table$V16 < 0] = 0
  histone.table$V16[is.na(histone.table$V16)] = 0
  histone.table.aggregate <- aggregate(list(readcov_1 = histone.table$V15, readcov_2 = histone.table$V16),
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
  colnames(count) = c(epi_id1, epi_id2)
  return(count)
}

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine2")
gene_id = data.frame(fread("gene_id", header=FALSE, col.names=c("gene_id")))
gene_id = as.data.table(str_split_fixed(gene_id$gene_id,':',2)[,1])
epi_id1 = "E003"
epi_id2 = "E004"
epi_id3 = "E005"
epi_id4 = "E094"
epi_id5 = "E095"
epi_id6 = "E096"

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

#===========================================
expression.count1.cor = get_cor_sample_pair(expression.count1, epi_id1, epi_id2)
expression.count2.cor = get_cor_sample_pair(expression.count2, epi_id1, epi_id3)
expression.count3.cor = get_cor_sample_pair(expression.count3, epi_id2, epi_id3)

methylation.count1.cor = get_cor_sample_pair(methylation.count1, epi_id1, epi_id2)
methylation.count2.cor = get_cor_sample_pair(methylation.count2, epi_id1, epi_id3)
methylation.count3.cor = get_cor_sample_pair(methylation.count3, epi_id2, epi_id3)


H3K36me3.count1.cor = get_cor_sample_pair(H3K36me3.count1, epi_id1, epi_id2)
H3K36me3.count2.cor = get_cor_sample_pair(H3K36me3.count2, epi_id1, epi_id3)
H3K36me3.count3.cor = get_cor_sample_pair(H3K36me3.count3, epi_id2, epi_id3)

#===========================================
expression.count1.cor.summary = summarize_cor(expression.count1.cor, 0.05, 0.7)
expression.count2.cor.summary = summarize_cor(expression.count2.cor, 0.05, 0.7)
expression.count3.cor.summary = summarize_cor(expression.count3.cor, 0.05, 0.7)

expression.count1.cor.summary = summarize_cor(expression.count1.cor, 1, 0)
expression.count2.cor.summary = summarize_cor(expression.count2.cor, 1, 0.7)
expression.count3.cor.summary = summarize_cor(expression.count3.cor, 1, 0.5)
lapply(expression.count2.cor.summary, function(x) length(x))

#================== NEW CORRELATION =====================
get_met_diff_cor <- function(epi_id1, epi_id2){
  print("---> Loading data...")
  file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  met.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  met.diff.exon = fread(met.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  met.diff.exon$V3 = "exon"
  met.diff.pi = fread(met.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_tables_list = list(met.diff.exon, met.diff.pi)
  
  for (i in 1:2){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table <- metdiff_table %>%
      mutate(V16 = as.numeric(as.character(V16)), V16 = if_else(is.na(V16), 0, V16)) %>%
      dplyr::select(V1, V3, V4, V5, V7, V9, V16)
    
    print("---> Changing id")
    if (i == 1){ #for exon table
      new_id = as.data.table(str_split_fixed(metdiff_table$V9, "\\;|\\ ", 8))[,c(2,8)]
      new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
      new_id = data.table(paste(new_id$V2,new_id$V8, sep=':E'))
      metdiff_table$V9 = new_id$V1
    }
    else {
      new_id = as.data.table(str_split_fixed(metdiff_table$V9, "\\; |\\ ", 4))[,c(2,4)]
      new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
      new_id = data.table(paste(new_id$V2,new_id$V4, sep=':I'))
      metdiff_table$V9 = new_id$V1
    }
    
    print("---> Aggregating...")
    metdiff_table.aggregate <- aggregate(list(metdiff = metdiff_table$V16), FUN=max,
                                   by = list(chr = metdiff_table$V1, feature = metdiff_table$V3,
                                             start = metdiff_table$V4, end = metdiff_table$V5, 
                                             strand = metdiff_table$V7, gene_id = metdiff_table$V9)
    )
    
    metdiff_tables_list[[i]] = metdiff_table.aggregate
  }
  exon_met = metdiff_tables_list[[1]]
  pi_met = metdiff_tables_list[[2]]
  # colnames(pi_met) = c("chr", "feature", "start","end", "strand", "id", "metdiff")

  metdiff = rbind(exon_met, pi_met)
  colnames(metdiff) = c("chr", "feature", "start","end", "strand", "id", "metdiff")
  metdiff = as(metdiff, "GRanges")
  metdiff = sort(metdiff)
  # return(metdiff_tables_list[[1]]) #Right now only compare exons
  return(metdiff)
}
pair1_2.metdiff = get_met_diff_cor(epi_id1, epi_id2)
pair1_3.metdiff = get_met_diff_cor(epi_id1, epi_id3)
pair2_3.metdiff = get_met_diff_cor(epi_id2, epi_id3)
pair4_5.metdiff = get_met_diff_cor(epi_id4, epi_id5)
pair4_6.metdiff = get_met_diff_cor(epi_id4, epi_id6)
pair5_6.metdiff = get_met_diff_cor(epi_id5, epi_id6)

get_his_diff_cor <- function(type, epi_id1, epi_id2){
  print("---> Loading data...")
  file_path = normalizePath(paste(getwd(), type, "annotatedcounts", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  his.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  his.diff.exon = fread(his.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  his.diff.exon$V3 = "exon"
  his.diff.pi = fread(his.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff_tables_list = list(his.diff.exon, his.diff.pi)
  for (i in 1:2){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table <- hisdiff_table %>%
      mutate(V13 = as.numeric(as.character(V13)), V13 = if_else(is.na(V13), 0, V13)) %>%
      dplyr::select(V1, V3, V4, V5, V7, V9, V13)
    print("---> Changing id")
    if (i == 1){ #for exon table
      new_id = as.data.table(str_split_fixed(hisdiff_table$V9, "\\;|\\ ", 8))[,c(2,8)]
      new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
      new_id = data.table(paste(new_id$V2,new_id$V8, sep=':E'))
      hisdiff_table$V9 = new_id$V1
    }
    else {
      new_id = as.data.table(str_split_fixed(hisdiff_table$V9, "\\; |\\ ", 4))[,c(2,4)]
      new_id = as.data.frame(sapply(new_id, function(x) gsub("\"", "", x)))
      new_id = data.table(paste(new_id$V2,new_id$V4, sep=':I'))
      hisdiff_table$V9 = new_id$V1
    }
    print("---> Aggregating...")
    hisdiff_table.aggregate <- aggregate(list(hisdiff = hisdiff_table$V13), FUN=max,
                                         by = list(chr = hisdiff_table$V1, feature = hisdiff_table$V3,
                                                   start = hisdiff_table$V4, end = hisdiff_table$V5, 
                                                   strand = hisdiff_table$V7, gene_id = hisdiff_table$V9)
    )
    hisdiff_tables_list[[i]] = hisdiff_table.aggregate
  }
  exon_his = hisdiff_tables_list[[1]]
  pi_his = hisdiff_tables_list[[2]]

  hisdiff = rbind(exon_his, pi_his)
  colnames(hisdiff) = c("chr", "feature", "start","end", "strand", "id", "hisdiff")
  hisdiff = as(hisdiff, "GRanges")
  hisdiff = sort(hisdiff)
  # return(hisdiff_tables_list[[1]]) #Right now only compare exons
  return(hisdiff)
}
pair1_2.hisdiff = get_his_diff_cor("H3K36me3", epi_id1, epi_id2)
pair1_3.hisdiff = get_his_diff_cor("H3K36me3", epi_id1, epi_id3)
pair2_3.hisdiff = get_his_diff_cor("H3K36me3", epi_id2, epi_id3)
pair4_5.hisdiff = get_his_diff_cor("H3K36me3", epi_id4, epi_id5)
pair4_6.hisdiff = get_his_diff_cor("H3K36me3", epi_id4, epi_id6)
pair5_6.hisdiff = get_his_diff_cor("H3K36me3", epi_id5, epi_id6)

get_exp_diff_cor <- function(epi_id1, epi_id2){
  print("---> Loading data...")
  file_path = normalizePath(paste(getwd(), "expression", "annotatedcounts", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  exp.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  exp.diff.exon = fread(exp.diff[[1]], header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
  exp.diff.exon = exp.diff.exon %>%
    mutate(gene_id = paste(groupID, featureID, sep=':')) %>%
    dplyr::select(c(8,7))
  colnames(exp.diff.exon) = c("gene_id", "log2fold")
  exp.diff.exon = exp.diff.exon %>% mutate(log2fold = if_else(is.na(log2fold), 0, log2fold))
  exp.diff.exon <- aggregate(list(log2fold = exp.diff.exon[[2]]), FUN=max,
                     by = list(gene_id = exp.diff.exon$gene_id))
  return(exp.diff.exon)
}
pair1_2.expdiff = get_exp_diff_cor(epi_id1, epi_id2)
pair1_3.expdiff = get_exp_diff_cor(epi_id1, epi_id3)
pair2_3.expdiff = get_exp_diff_cor(epi_id2, epi_id3)
pair4_5.expdiff = get_exp_diff_cor(epi_id4, epi_id5)
pair4_6.expdiff = get_exp_diff_cor(epi_id4, epi_id6)
pair5_6.expdiff = get_exp_diff_cor(epi_id5, epi_id6)

#======COR Exp and Methylation========
get_gene_cor <- function(expdiff, metdiff, hisdiff){
  print("---> Preparing exp_his, exp_met, met_his")
  exp_met = cbind(gene_id$V1,
                  expdiff$log2fold, 
                  as.data.frame(metdiff[metdiff$feature == "exon"]$metdiff))
  exp_his = cbind(gene_id$V1,
                  expdiff$log2fold,
                  as.data.frame(hisdiff[hisdiff$feature == "exon"]$hisdiff))

  met_his = cbind(str_split_fixed(metdiff$id, ':', 2)[,1],
                  as.data.frame(metdiff$metdiff),
                  as.data.frame(hisdiff$hisdiff))
  colnames(exp_met) = c("gene_id","log2fold", "diff")
  colnames(exp_his) = c("gene_id","log2fold", "diff")
  colnames(met_his) = c("gene_id","log2fold", "diff")

  diff_list = list(exp_met, exp_his, met_his)

  print("---> Calculating correlation for genes")
  cor_list = vector("list",3)
  for (i in 1:3){
    diff_table = diff_list[[i]]
    diff_table.cor = diff_table %>%
      group_by(gene_id) %>%
      mutate(cor = cor(log2fold, diff)) %>%
      dplyr::select(gene_id, cor) %>%
      unique()
    cor_list[[i]] = diff_table.cor
  }
  return(cor_list)
}
get_all_pair_cor <- function(pair_list){
  n_pair = length(pair_list)
  all_pair_cor = vector("list", n_pair)
  for (i in 1:n_pair){
    print("New pair")
    pair = pair_list[[i]]
    expdiff = pair[[1]]
    metdiff = pair[[2]]
    hisdiff = pair[[3]]
    all_pair_cor[[i]] = get_gene_cor(expdiff, metdiff, hisdiff)
  }
  return(all_pair_cor)
}

#for each pair, there are 3 groups: expmet, exphis, methis.
pair_list = list(list(pair1_2.expdiff, pair1_2.metdiff, pair1_2.hisdiff),
                 list(pair1_3.expdiff, pair1_3.metdiff, pair1_3.hisdiff),
                 list(pair2_3.expdiff, pair2_3.metdiff, pair2_3.hisdiff),
                 list(pair4_5.expdiff, pair4_5.metdiff, pair4_5.hisdiff),
                 list(pair4_6.expdiff, pair4_6.metdiff, pair4_6.hisdiff),
                 list(pair5_6.expdiff, pair5_6.metdiff, pair5_6.hisdiff))
all_pair_cor = get_all_pair_cor(pair_list) 

#Get dataframes of number of genes having cor/anticore between (1) exp_met (2) exp_his
all_pair_ncor.met = vector("list", 6) 
all_pair_ncor.his = vector("list", 6) 
all_pair_nanticor.met = vector("list", 6) 
all_pair_nanticor.his = vector("list", 6) 
for (i in 1:length(all_pair_cor)){
  pair = all_pair_cor[[i]]
  all_pair_ncor.met[[i]] = table(pair[[1]]$cor >= 0.5)[[2]]/length(unique(gene_id$V1))*100
  all_pair_ncor.his[[i]] = table(pair[[2]]$cor >= -0.5)[[2]]/length(unique(gene_id$V1))*100
  all_pair_nanticor.met[[i]] = table(pair[[1]]$cor <= -0.5)[[2]]/length(unique(gene_id$V1))*100
  all_pair_nanticor.his[[i]] = table(pair[[2]]$cor <= -0.5)[[2]]/length(unique(gene_id$V1))*100
  # print(count(pair[[1]]$cor <= -0.5)$freq)
  # print(count(pair[[1]]$cor > -0.5 & pair[[1]]$cor < 0.5)$freq)
}
all_pair_ncor.met = do.call(rbind.data.frame, all_pair_ncor.met)
all_pair_ncor.his = do.call(rbind.data.frame, all_pair_ncor.his)
all_pair_nanticor.met = do.call(rbind.data.frame, all_pair_nanticor.met)
all_pair_nanticor.his = do.call(rbind.data.frame, all_pair_nanticor.his)

#Get dataframes of genes having cor/anticore between (1) exp_met (2) exp_his
all_pair_corgenes = vector("list", 6)
for (i in 1:length(all_pair_cor)){
  pair = all_pair_cor[[i]]
  sig_genes = as.data.table(pair[[2]][pair[[2]]$cor > 0.7 & !is.na(pair[[2]]$cor), 1]$gene_id,
                            stringAsFactor = FALSE)
  # print(length(sig_genes$gene_id))
  all_pair_corgenes[[i]] = sig_genes
  print(sig_genes)
}

all_pair_corgenes = Reduce(union, all_pair_corgenes)
all_pair_corgenes = Reduce(union, all_pair_corgenes)
all_pair_corgenes = as.data.frame(all_pair_corgenes)
paste(all_pair_corgenes, collapse = ',')

known_genes = c("BRPF1", "DNMT3A", "GLYR1", "HDGF", "IWS1","MORF4L1"
                ,"MSH6","MTF2","MSL3","MUM1","NSD1","PHF1","PSIP1",
                "SPT16H","WHSC1","ZMYND11")
intersect(all_pair_corgenes, known_genes)




temp_cor = get_gene_cor(pair1_2.expdiff, pair1_2.hisdiff, pair1_2.metdiff)
pair1_2.exp_met.cor = temp_cor[[1]]
pair1_2.exp_his.cor = temp_cor[[2]]
pair1_2.met_his.cor = temp_cor[[3]]
length(pair1_2.exp_met.cor)

temp = as.data.table(pair1_2.exp_his.cor[pair1_2.exp_his.cor$cor > 0.7 & !is.na(pair1_2.exp_his.cor$cor), 1]$gene_id,
                     stringAsFactor = FALSE)
pair1_2.exp_his.cor = data.frame(pair1_2.exp_his.cor$gene_id, stringsAsFactors = FALSE)
temp2 = as.data.table(pair1_2.exp_his.cor[pair1_2.exp_his.cor$cor > 0.7 & !is.na(pair1_2.exp_his.cor$cor), 1]$gene_id,
                     stringAsFactor = FALSE)
temp = c(union(temp, temp2))
temp2 = as.data.table(pair1_2.exp_his.cor[pair1_2.exp_his.cor$cor > 0.7 & !is.na(pair1_2.exp_his.cor$cor),1])
temp = pair1_2.exp_his.cor %>%
  select(gene_id, cor) %>%
  filter(cor > 0.7) 
temp = as.data.frame(temp)
temp$gene_id



