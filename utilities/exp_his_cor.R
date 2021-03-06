library(RRreg)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
setwd("/home/dhthutrang/files/all_flank")
doMC::registerDoMC(cores = 17)
#===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
print("===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====")
all_pairs.exp = list.files("/home/dhthutrang/files/mRNA_seq/backupcount_NCBI/res", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp)$stat
  }
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  colnames(pair.exp_list) = c("gene_id", "exon_id", seq(1,(ncol(pair.exp_list)-2),1))
  return(pair.exp_list)
}
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ",i, sep=''))
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp)[,c("stat", "padj")]
  }
  pair.exp_list = lapply(pair.exp_list,
                         function(x) {
                           x = x %>%
                             mutate(exp = if_else(padj <= 0.05 & !is.na(padj), true = stat,
                                                  false = 0.0)) %>%
                             dplyr::select(exp)
                         })
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  colnames(pair.exp_list) = c("gene_id", "exon_id", seq(1,(ncol(pair.exp_list)-2),1))
  return(pair.exp_list)
}

#all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
#saveRDS(all_pairs.exp, "all_pairs.exp.RDS")
all_pairs.exp = readRDS("all_pairs.exp.RDS")
head(all_pairs.exp)
all_genes = fread("gene_id.txt", header = FALSE)

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
his_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    pair.his = pair.his %>%
      mutate(temp_val = abs(as.numeric(as.character(V10))),
             m_val = if_else(!is.na(temp_val) & temp_val >= 1, true = temp_val, 
                             false = 0)) %>%
      dplyr::select(m_val)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(max) %>%
    dplyr::select(-group)
  pair.his_list = cbind(his_id, pair.his_list)
  pair.his_list = pair.his_list[order(pair.his_list$V1)]
  return(pair.his_list)
}
get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    # for (j in 1:1){
    his = histone_type_list[[j]]
    print(his)
    all_pairs.his = list.files(his, full.names = TRUE)
    print(all_pairs.his)
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
    colnames(all_pairs.his.sig) = c("gene_id", "exon_id", seq(1, (ncol(all_pairs.his.sig)-2),1))
    all_pairs.his_list[[j]] = all_pairs.his.sig
  }
  return(all_pairs.his_list)
}
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac")
# all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
# head(all_pairs.his_list[[1]][[1]])
# saveRDS(all_pairs.his_list, "all_pairs.his_list.RDS")
all_pairs.his_list = readRDS("all_pairs.his_list.RDS")
head(all_pairs.his_list[[1]])

#===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
pearcor_p <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    p_val = p_value_calculator(cor(exp, his, method = "pearson"), nrow = length(exp))
    p_adj = p.adjust(p_val, method = "fdr", n=27247)
    return(p_adj)
  }
  else {
    return(NA)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.his, n_pairs){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  all_res_pair <- foreach( i=1:n_pairs, .combine='c', .packages=c('dplyr') ) %dopar% { #325 if other than H3K27ac and 231
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    res_table = data_table %>%
      group_by(all_pairs.exp$gene_id) %>%
      summarise(res = pearcor_p(exp, his)) %>%
      dplyr::select(res)
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method){
  all_res_list = vector("list", length(histone_type_list)-1 )
  for (j in 1:length(histone_type_list)){
    print(paste("Histone: ", histone_type_list[[j]], sep = ''))
    all_pairs.his = all_pairs.his_list[[j]]
    if (j == length(histone_type_list)){
      all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, n_pairs = 22*21/2)
    }
    else {
      all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, n_pairs = 25*24/2)
    }
    all_res_list[[j]] = all_res_pair
  }
  return(all_res_list)
}

print("Pearsons-p correlation")
all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, "pearcor_p")
saveRDS(all_res_list.pearcor_p, "all_res_list.pearcor_p.his.RDS")
