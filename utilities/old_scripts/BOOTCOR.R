library(RRreg)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
setwd("/home/dhthutrang/files/analyze_flank")
doMC::registerDoMC(cores = 17)
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
#===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
print("===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====")
all_pairs.exp = list.files("exp", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    pair.exp = all_pairs.exp[[i]]
    pair.exp = fread(pair.exp)[,6:7]
    colnames(pair.exp) = c("padj", "log2fold")
    pair.exp_list[[i]] = pair.exp
  }
  pair.exp_list = lapply(pair.exp_list,
                         function(x) {
                           colnames(x) = c("padj", "log2fold")
                           x = x %>%
                             mutate(exp = if_else(is.na(padj) | padj > 0.05, true = 0.0, 
                                                  false = if_else(log2fold > 0, true = 1.0,
                                                                  false = -1.0) ) ) %>%
                             dplyr::select(exp)
                         })
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  colnames(pair.exp_list) = c("gene_id", "exon_id", seq(1,(ncol(pair.exp_list)-2),1))
  return(pair.exp_list)
}

all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
# saveRDS(all_pairs.exp, "server_cat3/all_pairs.exp.RDS")
all_pairs.exp = readRDS("server_cat/all_pairs.exp.RDS")
head(all_pairs.exp)
all_genes = fread("gene_id.txt", header = FALSE)

#===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
  # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    pair.his = pair.his %>%
      mutate(m_val = as.numeric(as.character(V10)),
             p_val = as.numeric(as.character(V11)),
             sig = if_else( m_val >= 1 & p_val <= 0.05 & !is.na(p_val), true=1.0,
                            false = if_else(m_val <= -1 & p_val <= 0.05 & !is.na(p_val), true = -1,
                                            false = 0.0) )) %>%
      dplyr::select(sig)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(sum) %>%
    dplyr::select(-group)
  pair.his_list[pair.his_list < 0] = -1
  pair.his_list[pair.his_list > 0] = 1
  his_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
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
n_pairs = 15
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")
# all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
# saveRDS(all_pairs.his_list, "server_cat3/all_pairs.his_list.RDS")
all_pairs.his_list = readRDS("all_pairs.his_list.RDS")

#===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
randcor <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    randcor = RRcor(x = exp, y = his, models=c("Warner", "Warner"), p.list=list(0.2, 0.2))$r[[2]]
    return(randcor)
  }
  else {
    return(NA)
  }
}
pearcor <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(cor(exp, his, method = "pearson"))
  }
  else {
    return(NA)
  }
}
spearcor <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(cor(exp, his, method = "spearman"))
  }
  else {
    return(NA)
  }
}
pearcor_p <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(p_value_calculator(cor(exp, his, method = "pearson"), nrow = length(exp)))
  }
  else {
    return(NA)
  }
}
spearcor_p <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(p_value_calculator(cor(exp, his, method = "spearman"), nrow = length(exp)))
  }
  else {
    return(NA)
  }
}
fisher <- function(exp, his){
  exp = as.factor(exp)
  his = as.factor(his)
  if ( nlevels(exp) > 1 & nlevels(his) > 1){
    return(fisher.test(exp, his)$p.value)
  }
  else {
    return(NA)
  }
}
chisq <- function(exp, his){
  exp = as.factor(exp)
  his = as.factor(his)
  if ( nlevels(exp) > 1 & nlevels(his) > 1){
    return( chisq.test(exp, his)$p.value )
  }
  else {
    return(NA)
  }
}
theta_cor_fisher <- function(data,i){
  xdata <- data[i,]
  return(cor(as.numeric(xdata[, 1][[1]]), as.numeric(xdata[, 2][[1]])))
}
bootcor_fisher <- function(dtable, nperm, theta){
  if (length(unique(dtable[[1]])) > 1 & length(unique(dtable[[2]])) > 1 & nrow(dtable) > 3){
    # boot_cor <- bootstrap(x = seq_len(nrow(dtable)), nboot = nperm, theta, dtable)
    boot_cor <- boot(dtable, theta_cor_fisher, nperm, parallel = "multicore", ncpus = 4)$t
    mean_boot_cor = mean(boot_cor[!is.na(boot_cor)])
    # print(paste("New gene: ", mean_boot_cor, sep =''))
    if (!is.na(mean_boot_cor)){
      # print(paste("return mean_boot_cor: ", mean_boot_cor, sep =''))
      return(mean_boot_cor)
    }
    else{
      # print("return 2")
      return(2.0)
    }
  }
  else {
    # print("return 3")
    return(3.0)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.his, method, nperm = 1000, theta = theta_cor){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  # for (i in 1:(length(all_pairs.exp) - 2) ){
  all_res_pair <- foreach( i=1:(length(all_pairs.exp) - 2), .combine='c', .packages=c('dplyr') ) %dopar% {
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    
    if (method == "randcor"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = randcor(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "bootcor_fisher"){
      print(Sys.time())
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        do(
          summarise(., res = bootcor_fisher(dtable = ., nperm = nperm, theta = theta_cor_fisher))
        ) %>%
        ungroup %>%
        dplyr::select(res)
    }
    else if (method == "pearcor"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = pearcor(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "pearcor_p"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = pearcor_p(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "spearcor"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = spearcor(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "spearcor_p"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = spearcor_p(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "fisher"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = fisher(exp, his)) %>%
        dplyr::select(res)
    }
    else if (method == "chisq"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = chisq(exp, his)) %>%
        dplyr::select(res)
    }
    # all_res_pair[[i]] = res_table
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method){
  all_res_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
    print(paste("Histone: ", histone_type_list[[j]], sep = ''))
    all_pairs.his = all_pairs.his_list[[j]]
    all_res_pair = analyze_array(all_pairs.exp, all_pairs.his, method)
    all_res_list[[j]] = all_res_pair
  }
  return(all_res_list)
}
# 
# print("Pearsons correlation")
# all_res_list.pearcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "pearcor")
# saveRDS(all_res_list.pearcor, "all_res_list.pearcor.RDS")
# print("Spearman correlation")
# all_res_list.spearcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "spearcor")
# saveRDS(all_res_list.spearcor, "all_res_list.spearcor.RDS")
print("Pearsons_p correlation")
all_res_list.pearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, "pearcor_p")
saveRDS(all_res_list.pearcor_p, "all_res_list.pearcor_p.RDS")
print("Spearman_p correlation")
all_res_list.spearcor_p = analyze_array_list(all_pairs.exp, all_pairs.his_list, "spearcor_p")
saveRDS(all_res_list.spearcor_p, "all_res_list.spearcor_p.RDS")
# print("Randomized correlation")
# all_res_list.randcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "randcor")
# saveRDS(all_res_list.randcor, "all_res_list.randcor.RDS")
# print("Fisher Exact Test")
# all_res_list.fisher = analyze_array_list(all_pairs.exp, all_pairs.his_list, "fisher")
# saveRDS(all_res_list.fisher, "all_res_list.fisher.RDS")
# print("Chi Square")
# all_res_list.chisq = analyze_array_list(all_pairs.exp, all_pairs.his_list, "chisq")
# saveRDS(all_res_list.chisq, "all_res_list.chisq.RDS")
# print("Bootstrapped correlation")
# all_res_list.bootcor_fisher = analyze_array_list(all_pairs.exp, all_pairs.his_list, "bootcor_fisher")
# saveRDS(all_res_list.bootcor_fisher, "all_res_list.bootcor_fisher.RDS")
# print(gsub(":", "-", Sys.time()))
