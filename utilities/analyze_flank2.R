library(RRreg)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(bootstrap)
library(boot)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")

#===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
all_pairs.exp = list.files("exp", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp)$padj
  }
  pair.exp_list = lapply(pair.exp_list, 
                         function(x) ifelse((x > 0 & x <= 0.05 & !is.na(x)), 1.0, 0.0))
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  colnames(pair.exp_list) = c("gene_id", "exon_id", seq(1,(ncol(pair.exp_list)-2),1))
  return(pair.exp_list)
}

all_pairs.exp = get_all_pairs.exp(all_pairs.exp)
head(all_pairs.exp)
all_genes = fread("gene_id.txt", header = FALSE)


#===== PREPARE HIS FILE (6 TOTAL) =====
get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    pair.his = pair.his %>%
      mutate(m_val = as.numeric(as.character(V10)),
             p_val = as.numeric(as.character(V11)),
             sig = if_else( abs(m_val) >= 1 & p_val <= 0.05 & !is.na(p_val), 1.0, 0.0)) %>%
      select(sig)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(min) %>%
    select(-group)
  his_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.his_list = cbind(his_id, pair.his_list)
  pair.his_list = pair.his_list[order(pair.his_list$V1)]
  return(pair.his_list)
}
get_all_pairs.his_list <- function(histone_type_list){
  all_pairs.his_list = vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)){
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
all_pairs.his_list = get_all_pairs.his_list(histone_type_list)


#===== CORRELATION WITH RANDOMIZATION =====
temp.exp = all_pairs.exp$`1`
temp.his = all_pairs.his[[1]]$`1`
temp.cor = RRcor(x = temp.exp, y = temp.his, models=c("Warner", "Warner"), p.list=list(0.2, 0.2))
temp.cor$r[[2]]
cor(temp.exp, temp.his)
fisher.test(temp.exp, temp.his)$p
temp = sample(temp.exp)
table(temp == temp.exp)
cor(temp.exp, temp.his)
cor(temp, temp.his)

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
    return(cor(exp, his))
  }
  else {
    return(NA)
  }
}
fisher <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(fisher.test(exp, his)$p)
  }
  else {
    return(NA)
  }
}

theta_cor <- function(x, xdata){ 
  p_value_calculator(cor(as.numeric(xdata[x, 2][[1]]), as.numeric(xdata[x, 3][[1]])), nrow(xdata))
}

bootcor <- function(dtable, nperm, theta){
  cor_r = cor(as.numeric(dtable[[1]]), as.numeric(dtable[[2]]))
  # print(cor_r)
  if ( abs(cor_r) >= 0.5 & !is.na(cor_r) ){
    boot_cor <- bootstrap(x = seq_len(nrow(dtable)), nboot = nperm, theta, dtable)
    boot_cor = boot_cor$thetastar
    bootp = nperm*0.05
    # print(bootp)
    # print(length(boot_cor[boot_cor <= 0.05 & !is.na(boot_cor)]))
    # print(boot_cor[boot_cor <= 0.05 & !is.na(boot_cor)])
    if (length(boot_cor[boot_cor <= 0.05 & !is.na(boot_cor)]) <= bootp){
      return(1.0)
    } #1: less than 5% of all permutation has significant p-value
    else{
      return(0.0)
    }
  }
}
analyze_array <- function(all_pairs.exp, all_pairs.his, method, nperm = 100, theta = theta_cor){
  all_res_pair = vector("list", ncol(all_pairs.exp) - 2)
  for (i in 1:(length(all_pairs.exp) - 2) ){
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp[[i+2]]
    his = all_pairs.his[[i+2]]
    data_table = as.data.table(cbind(exp, his))
    if (method == "randcor"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = randcor(exp, his)) %>%
        select(res)
    }
    else if (method == "bootcor"){
      res_table = data_table %>%
        dplyr::group_split(all_pairs.exp$gene_id)
      res_table_boot = lapply(res_table, function(x) bootcor(dtable = x, nperm = 100, theta = theta_cor))
      names(res_table_boot) = all_genes$V1
      # res_table = t(as.data.table(res_table_boot))
      res_table = colnames(as.data.table(res_table_boot))
    }
    else if (method == "pearcor"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = pearcor(exp, his)) %>%
        select(res)
    }
    else if (method == "fisher"){
      res_table = data_table %>%
        group_by(all_pairs.exp$gene_id) %>%
        summarise(res = fisher(exp, his)) %>%
        select(res)
    }
    all_res_pair[[i]] = res_table
  }
  if (method != "bootcor"){
    all_res_pair = as.data.table(all_res_pair)
    all_res_pair = cbind(all_genes, all_res_pair)
    return(all_res_pair)
  }
  else {
    return(all_res_pair)
  }
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
# all_res_list.pearcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "pearcor")
# all_res_list.randcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "randcor")
# all_res_list.fisher = analyze_array_list(all_pairs.exp, all_pairs.his_list, "fisher")
# all_res_list.bootcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "bootcor")

# ------------ Get significant results ------------
get_all_res_list_sig <- function(all_res_list, method){
  all_res_list_sig = vector("list", length(all_res_list))
  for (i in 1:length(all_res_list)) {
    all_res = all_res_list[[i]]
    all_res_sig = vector("list", ncol(all_res)-1)
    for (j in 2:ncol(all_res)){
      if (method == "randcor"){
        all_res_sig[j-1] = all_res[abs(all_res[[j]]) >= 0.7  & !is.na(all_res[[j]]), 1] 
      }
      else if (method == "pearcor"){
        all_res_sig[j-1] = all_res[abs(all_res[[j]]) >= 0.7  & !is.na(all_res[[j]]), 1] 
      }
      else if (method == "fisher"){
        all_res_sig[j-1] = all_res[all_res[[j]] <= 0.05  & !is.na(all_res[[j]]), 1]
      }
    }
    all_res_list_sig[[i]] = all_res_sig
  }
  return(all_res_list_sig)
}
# all_res_list.pearcor_sig = get_all_res_list_sig(all_res_list.pearcor, "pearcor")
# all_res_list.randcor_sig = get_all_res_list_sig(all_res_list.randcor, "randcor")
# all_res_list.fisher_sig = get_all_res_list_sig(all_res_list.fisher, "fisher")
# all_res_list.bootcor_sig = all_res_list.bootcor

join_all_res_list_sig <- function(all_res_list_sig){
  all_res_list_sig_joined = vector("list")
  for (j in 1:n_pairs){
    all_res_sig = vector("list")
    for (i in 1:length(histone_type_list)){
      all_res_sig[[i]] = all_res_list_sig[[i]][[j]]
    }
    all_res_sig = Reduce(union, all_res_sig)
    print(length(all_res_sig))
    all_res_list_sig_joined[[j]] = all_res_sig
  }
  return(all_res_list_sig_joined)
}
# all_res_list.pearcor_sig_joined = join_all_res_list_sig(all_res_list.pearcor_sig)
# all_res_list.randcor_sig_joined = join_all_res_list_sig(all_res_list.randcor_sig)
# all_res_list.fisher_sig_joined = join_all_res_list_sig(all_res_list.fisher_sig)
# all_res_list.bootcor_sig_joined = join_all_res_list_sig(all_res_list.bootcor_sig)

get_tissue_specific_gene <- function(all_epigenes_list, pair_list){
  all_epigenes_intersect = Reduce(intersect, all_epigenes_list[pair_list])
  all_epigenes_intersect = gsub("\\+", ", ", all_epigenes_intersect)
  all_epigenes_intersect = paste(all_epigenes_intersect, collapse = ", ")
  all_epigenes_intersect = strsplit(all_epigenes_intersect, split=', ')[[1]]
  return(all_epigenes_intersect)
}
all_tissues = c(1:15) #1
esc = c(1:5) #Esc 2
mes = c(1,6,7,8,9) #Mesendoderm 3
tro = c(2,6,10,11,12) #Trophoblast 4 
gas = c(3,7,10,13,14) #Gastric 5 
ven = c(4,8,11,13,15) #Ventricle 6
mus = c(5,9,12,14,15) #Muscles 7
tissue_type_list = list(all_tissues, esc, mes, tro, gas, ven, mus)

all_res_list.pearcor_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.pearcor_sig_joined, x))
all_res_list.randcor_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.randcor_sig_joined, x))
all_res_list.fisher_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.fisher_sig_joined, x))
all_res_list.bootcor_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.bootcor_sig_joined, x))
paste(all_res_list.bootcor_sig_joined_genes[[5]], collapse = ', ')


##===== COMBINE =====
