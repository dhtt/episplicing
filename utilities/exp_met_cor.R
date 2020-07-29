library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(boot)
library(stats)
library(parallel)
library("doMC")
setwd("/home/dhthutrang/files/all_flank")
doMC::registerDoMC(cores = 17)
#===== PREPARE EXP FILE (1 FOR ALL met TYPES) =====
print("===== PREPARE EXP FILE (1 FOR ALL met TYPES) =====")
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
#===== PREPARE met FILE (6 TOTAL) =====
met_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
get_all_pairs.met <- function(all_pairs.met){
  pair.met_list = vector("list", length(all_pairs.met))
  for (i in 1:length(all_pairs.met)){
    print(paste("Pair: ", i, sep=''))
    pair.met = all_pairs.met[[i]]
    pair.met = fread(pair.met)
    pair.met = pair.met %>%
      mutate(temp_val = abs(as.numeric(as.character(V10))),
             m_val = if_else(!is.na(temp_val) & temp_val >= 25, true = temp_val,
                             false = 0)) %>%
      dplyr::select(m_val)
    pair.met_list[[i]] = pair.met
  }
  pair.met_list = as.data.table(pair.met_list)
  print("======== GROUPING ======== ")
  pair.met_list = pair.met_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(max) %>%
    dplyr::select(-group)
  pair.met_list = cbind(met_id, pair.met_list)
  print("======== ORDERING ======== ")
  pair.met_list = pair.met_list[order(pair.met_list$V1)]
  return(pair.met_list)
}

get_all_pairs.met_list <- function(){
  met = "met"
  print(met)
  all_pairs.met = list.files(met, full.names = TRUE)
  print(all_pairs.met)
  all_pairs.met.sig = get_all_pairs.met(all_pairs.met)
  colnames(all_pairs.met.sig) = c("gene_id", "exon_id", seq(1, (ncol(all_pairs.met.sig)-2),1))
  return(all_pairs.met.sig)
}
#all_pairs.met.sig = get_all_pairs.met_list()
#head(all_pairs.met.sig)
#saveRDS(all_pairs.met.sig, "all_pairs.met_list.RDS")
all_pairs.met_list = readRDS("all_pairs.met_list.RDS")
head(all_pairs.met_list)

#===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
pearcor_p <- function(exp, met){
  if (length(unique(exp)) > 1 & length(unique(met)) > 1){
    p_val = p_value_calculator(cor(exp, met, method = "pearson"), nrow = length(exp))
    p_adj = p.adjust(p_val, method = "fdr", n=27247)
    return(p_adj)
  }
  else {
    return(NA)
  }
}

analyze_array <- function(all_pairs.exp, all_pairs.met, n_pairs){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  all_res_pair <- foreach( i=1:n_pairs, .combine='c', .packages=c('dplyr') ) %dopar% { #325 if other than H3K27ac and 231
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp[[i+2]]
    met = all_pairs.met[[i+2]]
    data_table = as.data.table(cbind(exp, met))
    res_table = data_table %>%
      group_by(all_pairs.exp$gene_id) %>%
      summarise(res = pearcor_p(exp, met)) %>%
      dplyr::select(res)
  }
  all_res_pair = as.data.table(all_res_pair)
  all_res_pair = cbind(all_genes, all_res_pair)
  return(as.data.frame(all_res_pair))
}

print("Pearsons-p correlation")
all_res_list.pearcor_p = analyze_array(all_pairs.exp = all_pairs.exp, all_pairs.met = all_pairs.met_list, n_pairs = 300)
saveRDS(all_res_list.pearcor_p, "all_res_list.pearcor_p.met.RDS")

#temp = readRDS("/Users/dhthutrang/all_pairs.met_list.RDS")
#temp1 = readRDS("/Users/dhthutrang/all_pairs.his_list.RDS")
#temp1 = temp1[[1]]
#head(temp1[, 1:10])

