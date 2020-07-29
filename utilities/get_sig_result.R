library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(stringr)
library(boot)
library(stats)
library(parallel)
library("doMC")
setwd("/home/dhthutrang/files/all_flank")
doMC::registerDoMC(cores = 17)

print("===== PREPARE EXP FILE =====")
all_pairs.exp = list.files("/home/dhthutrang/files/mRNA_seq/backupcount_NCBI/res", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ",i, sep=''))
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp, header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
  }
  all_res_pair <- foreach( i=1:length(pair.exp_list), .combine='append', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    exp_table = pair.exp_list[[i]]
    exp_res = exp_table %>%
      filter(padj <= 0.05) %>%
      group_by(groupID) %>%
      mutate(DEU = paste(featureID, collapse = ', ')) %>%
      dplyr::select(groupID, DEU) %>%
      unique()
  }
  return(all_res_pair)
}
# all_DEU_res = get_all_pairs.exp(all_pairs.exp)
# saveRDS(all_DEU_res, "all_DEU_res.exp.RDS")

print("===== PREPARE HIS FILE =====")
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    print(paste("Pair: ",i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his_list[[i]] = fread(pair.his, header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  }
  all_res_pair <- foreach( i=1:length(pair.his_list), .combine='append', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    his_table = pair.his_list[[i]]
    
    his_res = his_table %>%
      mutate(m_val = as.numeric(as.character(V10))) %>%
      filter(abs(m_val) >= 1) %>%
      group_by(V9) %>%
      dplyr::select(V9) %>%
      unique()
    his_res = as.data.table(str_split_fixed(his_res$V9, pattern = '"', 8)[,c(2,6)])
    his_res = his_res %>%
      group_by(V1) %>%
      mutate(featureID = paste(V2, collapse = ', ')) %>%
      dplyr::select(V1, featureID) %>%
      unique()
  }
  return(all_res_pair)
}
get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    his = histone_type_list[[j]]
    print(his)
    all_pairs.his = list.files(his, full.names = TRUE)
    print(all_pairs.his)
    all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
    all_pairs.his_list[[j]] = all_pairs.his.sig
  }
  return(all_pairs.his_list)
}

histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac")
# all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
# saveRDS(all_pairs.his_list, "all_DEU_res.his.RDS")

print("===== PREPARE MET FILE =====")
all_pairs.met = list.files("/home/dhthutrang/files/all_flank/met", full.names = TRUE)
get_all_pairs.met <- function(all_pairs.met){
  pair.met_list = vector("list", length(all_pairs.met))
  for (i in 1:length(all_pairs.met)){
    print(paste("Pair: ",i, sep=''))
    pair.met = all_pairs.met[[i]]
    pair.met_list[[i]] = fread(pair.met, header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  }
  all_res_pair <- foreach( i=1:length(pair.met_list), .combine='append', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    met_table = pair.met_list[[i]]
    
    met_res = met_table %>%
      filter(abs(V10) >= 25) %>%
      group_by(V9) %>%
      dplyr::select(V9) %>%
      unique()
    met_res = as.data.table(str_split_fixed(met_res$V9, pattern = '"', 8)[,c(2,6)])
    met_res = met_res %>%
      group_by(V1) %>%
      mutate(featureID = paste(V2, collapse = ', ')) %>%
      dplyr::select(V1, featureID) %>%
      unique()
  }
  return(all_res_pair)
}
# 
# all_DEU_res.met = get_all_pairs.met(all_pairs.met)
# saveRDS(all_DEU_res.met, "all_DEU_res.met.RDS")






























