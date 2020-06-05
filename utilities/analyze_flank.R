library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(Hmisc)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")

p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
get_tissue_specific_gene <- function(all_epigenes_list, pair_list){
  all_epigenes_intersect = Reduce(intersect, all_epigenes_list[pair_list])
  all_epigenes_intersect = gsub("\\+", ", ", all_epigenes_intersect)
  all_epigenes_intersect = paste(all_epigenes_intersect, collapse = ", ")
  all_epigenes_intersect = strsplit(all_epigenes_intersect, split=', ')[[1]]
  return(all_epigenes_intersect)
}
join_all_his_res <- function(all_epigenes_chi_list){
  all_res_list = vector("list")
  for (j in 1:15){
    all_res = vector("list")
    for (i in 1:length(histone_type_list)){
      all_res[[i]] = all_epigenes_chi_list[[i]][[j]]
    }
    all_res = Reduce(union, all_res)
    print(length(all_res))
    all_res_list[[j]] = all_res
  }
  return(all_res_list)
}
all_tissues = c(1:15)
esc = c(1:5) #Esc
mes = c(1,6,7,8,9) #Mesendoderm
tro = c(2,6,10,11,12) #Trophoblast
gas = c(3,7,10,13,14) #Gastric
ven = c(4,8,11,13,15) #Ventricle
mus = c(5,9,12,14,15) #Muscles
tissue_type_list = list(all_tissues, esc, mes, tro, gas, ven, mus)

#===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
all_pairs.exp = list.files("exp", full.names = TRUE)
get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp)$padj
  }
  pair.exp_list = lapply(pair.exp_list, function(x) ifelse((is.na(x)), 1.0, x))
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  return(pair.exp_list)
}
all_pairs.exp.padj = get_all_pairs.exp(all_pairs.exp)
colnames(all_pairs.exp.padj) = c("gene_id", "exon_id", seq(1,(ncol(all_pairs.exp.padj)-2),1))
head(all_pairs.exp.padj)
all_genes = fread("gene_id.txt", header = FALSE)

# flank = import.gff("reference_genome.fl300.gtf")
# flank_id = unique(as.data.table(cbind(flank$gene, flank$exon)))
# fwrite(flank_id, "flank_id.txt", sep = '\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

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
             sig = if_else( abs(m_val) >= 1 & p_val <= 0.05 & !is.na(p_val), p_val, 1.0)) %>%
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

histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")
all_pairs.his.sig_list = vector("list", length(histone_type_list))
for (i in 1:length(histone_type_list)){
  his = histone_type_list[[i]]
  print(his)
  all_pairs.his = list.files(his, full.names = TRUE)
  print(all_pairs.his)
  all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
  colnames(all_pairs.his.sig) = c("gene_id", "exon_id", seq(1, (ncol(all_pairs.his.sig)-2),1))
  all_pairs.his.sig_list[[i]] = all_pairs.his.sig
}

#============================Correlation============================
new_cor <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(cor(exp, his))
  }
  else {
    return(0.0)
  }
}
get_correlation <- function(all_pairs.his.sig, method){
  all_cors_list = vector("list", ncol(all_pairs.exp.padj)-2)
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp.padj[[i+2]]
    his = all_pairs.his.sig[[i+2]]
    cor_table = as.data.table(cbind(exp, his))
    cor_table = cor_table %>%
      group_by(all_pairs.exp.padj$gene_id) %>%
      summarise(p_adj = p.adjust(p_value_calculator(cor(exp, his, method = method), n())), "fdr", n()) %>%
                # r_val = cor(exp, his, method = method),
                # p_val = p_value_calculator(r_val, n()),
                # p_adj = p.adjust(p_val, "fdr", n())) 
      select(p_adj)
    all_cors_list[[i]] = cor_table
  }
  all_cors_list = as.data.table(all_cors_list)
  all_cors_list = cbind(all_genes, all_cors_list)
  return(all_cors_list)
}

all_cors_pearson_list = vector("list", length(all_pairs.his.sig_list))
for (i in 1:length(all_pairs.his.sig_list)){
  print(paste("Histone: ", histone_type_list[[i]], sep = ''))
  all_pairs.his.sig = all_pairs.his.sig_list[[i]]
  all_cors_pearson = get_correlation(all_pairs.his.sig, "pearson")
  all_cors_pearson_list[[i]] = all_cors_pearson
}

# all_p_val = as.data.frame(all_cors)[, c(1, seq(0, ncol(all_cors), 3))]
# all_p_adj_pearson = as.data.frame(all_cors_pearson)[, seq(1, ncol(all_cors_pearson), 3)]
# all_r_val = as.data.frame(all_cors)[, c(1, seq(2, ncol(all_cors), 3))]

all_epigenes_cor_pearson_list = vector("list", length(all_cors_pearson_list))
for (i in 1:length(all_cors_pearson_list)) {
  all_cors_pearson = all_cors_pearson_list[[i]]
  all_epigenes_cor_pearson = vector("list", ncol(all_cors_pearson)-1)
  for (j in 2:ncol(all_cors_pearson)){
    all_epigenes_cor_pearson[j-1] = all_cors_pearson[all_cors_pearson[[j]] <= 0.05  & !is.na(all_cors_pearson[[j]]), 1] 
  }
  all_epigenes_cor_pearson_list[[i]] = all_epigenes_cor_pearson
}

all_epigenes_cor_pearson_list = join_all_his_res(all_epigenes_cor_pearson_list)
all_epigenes_cor_pearson_list = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_cor_pearson_list, x))
all_epigenes_cor_pearson_list = transpose(as.data.table(all_epigenes_cor_pearson_list))
saveRDS(all_epigenes_cor_pearson_list, "all_epigenes_cor_pearson_list.RDS")
fwrite(temp, "cor_pearson_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')

#============================Correlation=========================
all_cors_spearman_list = vector("list", length(all_pairs.his.sig_list))
for (i in 1:length(all_pairs.his.sig_list)){
  print(paste("Histone: ", histone_type_list[[i]], sep = ''))
  all_pairs.his.sig = all_pairs.his.sig_list[[i]]
  all_cors_spearman = get_correlation(all_pairs.his.sig, "spearman")
  all_cors_spearman_list[[i]] = all_cors_spearman
}

# all_p_val = as.data.frame(all_cors)[, c(1, seq(0, ncol(all_cors), 3))]
# all_p_adj_spearman = as.data.frame(all_cors_spearman)[, seq(1, ncol(all_cors_spearman), 3)]
# all_r_val = as.data.frame(all_cors)[, c(1, seq(2, ncol(all_cors), 3))]

all_epigenes_cor_spearman_list = vector("list", length(all_cors_spearman_list))
for (i in 1:length(all_cors_spearman_list)) {
  all_cors_spearman = all_cors_spearman_list[[i]]
  all_epigenes_cor_spearman = vector("list", ncol(all_cors_spearman)-1)
  for (j in 2:ncol(all_cors_spearman)){
    all_epigenes_cor_spearman[j-1] = all_cors_spearman[all_cors_spearman[[j]] <= 0.05  & !is.na(all_cors_spearman[[j]]), 1] 
  }
  all_epigenes_cor_spearman_list[[i]] = all_epigenes_cor_spearman
}

all_epigenes_cor_spearman_list = join_all_his_res(all_epigenes_cor_spearman_list)
all_epigenes_cor_spearman_list = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_cor_spearman_list, x))
saveRDS(all_epigenes_cor_spearman_list, "all_epigenes_cor_spearman_list.RDS")
fwrite(all_epigenes_cor_spearman_list, "cor_spearman_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')
#============================Get Chi============================
new_chi.test <- function(exp, his){
  exp = as.factor(exp)
  his = as.factor(his)
  if ( nlevels(exp) > 1 & nlevels(his) > 1){
    return( chisq.test(exp, his)$p.value )
  }
  else {
    return(1.0)
  }
}
get_chi <- function(all_pairs.his.sig){
  all_chis_list = vector("list", ncol(all_pairs.exp.padj)-2)
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp.padj[[i+2]]
    his = all_pairs.his.sig[[i+2]]
    chi_table = as.data.table(cbind(exp, his))
    chi_table = chi_table %>%
      group_by(all_pairs.exp.padj$gene_id) %>%
      summarise(p_val = new_chi.test(exp, his)) %>%
      dplyr::select(p_val)
    all_chis_list[[i]] = chi_table
  }
  all_chis_list = as.data.table(all_chis_list)
  all_chis_list = cbind(all_genes, all_chis_list)
  return(all_chis_list)
}

all_chis_list = vector("list", length(all_pairs.his.sig_list))
for (i in 1:length(all_pairs.his.sig_list)){
  print(paste("Histone: ", histone_type_list[[i]], sep = ''))
  all_pairs.his.sig = all_pairs.his.sig_list[[i]]
  all_chis = get_chi(all_pairs.his.sig)
  all_chis_list[[i]] = all_chis
}

all_epigenes_chi_list = vector("list", length(all_chis_list))
for (i in 1:length(all_chis_list)) {
  all_chis = all_chis_list[[i]]
  all_epigenes_chi = vector("list", ncol(all_chis)-1)
  for (j in 2:ncol(all_chis)){
    all_epigenes_chi[j-1] = all_chis[all_chis[[j]] <= 0.05  & !is.na(all_chis[[j]]), 1] 
  }
  all_epigenes_chi_list[[i]] = all_epigenes_chi
}

temp = Reduce(union, all_epigenes_chi_list)
temp = Reduce(union, temp)
all_chis = temp
# saveRDS(temp, "all_epigenes_chi_list_unsplit.RDS")

all_epigenes_chi_list = join_all_his_res(all_epigenes_chi_list)
all_epigenes_chi_list = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_chi_list, x))
# saveRDS(all_epigenes_chi_list, "all_epigenes_chi_list.RDS")
fwrite(all_epigenes_chi_list, "chi_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')
#===========================EXCLUSIVE EXON===========================
get_all_exon_freq <- function(){
  all_sums_list = vector("list", length(all_pairs.his.sig_list))
  for (i in 1: length(all_pairs.his.sig_list)){
    all_pairs.his.sig = all_pairs.his.sig_list[[i]]
    combine_exp_his = as.data.frame(all_pairs.exp.padj <= 0.05 & all_pairs.his.sig < 1)[,3:17]
    
    all_sums = vector("list", length(tissue_type_list))
    for (j in 1:length(tissue_type_list)){
      tissue = tissue_type_list[[j]]
      sums = combine_exp_his[tissue]
      all_sums[[j]] = rowSums(sums)
    }
    all_sums = as.data.table(all_sums) 
    colnames(all_sums) = c("all", "esc", "mes", "tro", "gas", "ven", "mus")
    all_sums$gene_id = all_pairs.exp.padj$gene_id
    all_sums$exon_id = all_pairs.exp.padj$exon_id
    all_sums_list[[i]] = all_sums
  }
  return(all_sums_list)
}
all_exon_freqs = get_all_exon_freq()

get_all_sig_exon <- function(all_exon_freqs, threshold = 0){
  all_exons_list = vector("list", length(histone_type_list))
  for (i in 1:length(histone_type_list)){
    histone_type = all_exon_freqs[[i]]
    all_exons = vector("list", length(tissue_type_list))
    for (j in 1:length(tissue_type_list)){
      temp1 = unique(histone_type[histone_type[[j]] >= (max(histone_type[[j]] - threshold)), `gene_id`])
      all_exons[[j]] = temp1
    }
    all_exons_list[[i]] = all_exons
  }
  return(all_exons_list)
}
all_exon_list = get_all_sig_exon(all_exon_freqs, threshold = 0)

join_all_exon_res <- function(all_exons_list){
  all_res_list = vector("list")
  for (j in 1:length(tissue_type_list)){
    all_res = vector("list")
    for (i in 1:length(histone_type_list)){
      all_res[[i]] = all_exons_list[[i]][[j]]
    }
    all_res = Reduce(union, all_res)
    print(length(all_res))
    all_res_list[[j]] = all_res
  }
  return(all_res_list)
}
all_exon_list = join_all_exon_res(all_exon_list)

get_tissue_specific_exon <- function(all_exons){
  all_exons = gsub("\\+", ", ", all_exons)
  all_exons = paste(all_exons, collapse = ", ")
  all_exons = strsplit(all_exons, split=', ')[[1]]
  return(all_exons)
}
all_exon_list = lapply(all_exon_list, function(x) get_tissue_specific_exon(x))
# saveRDS(all_exon_list, "all_exon_list.RDS")
all_exon_list[[1]]
all_epigenes_chi_list[[1]]

library("VennDiagram", quietly=TRUE)
for (i in 1:length(all_exon_list)){
  venn.plot <- venn.diagram(x = list(
    "Gene-level (Chi-square)" = all_epigenes_chi_list[[i]],
    "Exon-level " = all_exon_list[[i]]
    # "From Chi-square" = all_exon_list[[1]]
  ), 
  filename=NULL, col="black", fill = c("yellow", "red"), margin=0.1, alpha = 0.5,
  cat.dist = 0.1, print.mode="percent")
  jpeg(paste(i, "_gene_exon.jpg"), res = 300, width = 7, height = 5, units="in"); 
  grid.draw(venn.plot); 
  dev.off()
}
#===========================EXCLUSIVE EXON COR===========================
get_all_exon_cor <- function(){
  all_sums_list = vector("list", length(all_pairs.his.sig_list))
  for (i in 1: length(all_pairs.his.sig_list)){
  # for (i in 1:1){
    print(paste("His: ", i, sep=''))
    all_pairs.his.sig = all_pairs.his.sig_list[[i]]
    combine_exp_his = as.data.table(all_pairs.exp.padj <= 0.05 & all_pairs.his.sig < 1)
    combine_exp_his[combine_exp_his == TRUE] = 1.0
    combine_exp_his[combine_exp_his == FALSE] = 0.0 
    combine_exp_his$gene_id = all_pairs.exp.padj$gene_id
    combine_exp_his$exon_id = all_pairs.exp.padj$exon_id
    combine_exp_his = combine_exp_his %>%
      dplyr::group_split(gene_id)
    combine_exp_his = lapply(combine_exp_his, function(x) exon_cor(x))
    combine_exp_his = unlist(combine_exp_his, use.names=FALSE)
    all_sums_list[[i]] = combine_exp_his
  }
  return(all_sums_list)
}
exon_cor <- function(all_exons){
  all_exons = all_exons %>%
    dplyr::select(-c(gene_id, exon_id)) %>%
    t %>%
    as.matrix %>%
    cor
  if ( length(table(abs(all_exons) >= 0.5 & abs(all_exons) < 1 & !is.na(all_exons))) > 1)
    return(1)
  else
    return(0)
}
exon_cor <- function(all_exons){
  all_exons = t(as.matrix(all_exons[,3:17]))
  all_exons = p_value_calculator(cor(all_exons), ncol(all_exons))
  if ( length(table(all_exons <= 0.05 & all_exons > 0 & !is.na(all_exons))) > 1)
    return(1)
  else
    return(0)
}
all_exon_cors = get_all_exon_cor()
# saveRDS(all_exon_cors, "all_exon_cors.RDS")
temp = all_exon_cors
temp = lapply(temp, as.data.table)
temp = lapply(temp, function(x) cbind(all_genes, x))
all_exons_res = vector("list", length(temp))
for (i in 1:length(temp)){
  cor_genes = temp[[i]]
  cor_genes = cor_genes[cor_genes[[2]] == 1, 1]
  # cor_genes = lapply(cor_genes, function(x) paste(x, collapse = ', '))
  all_exons_res[[i]] = cor_genes
}
all_exons_res = Reduce(union, all_exons_res)
temp = lapply(all_exons_res, function(x) gsub("\\+", ", ", x))
temp = lapply(temp, function(x)  paste(x, collapse = ', '))[[1]]
temp = strsplit(temp, split = ', ')[[1]]
exon_res_p = temp
exon_res = readRDS("exon_res_cor.RDS")
saveRDS(exon_res, "exon_res_cor.RDS")

#===========================EXCLUSIVE EXON FISHER===========================
get_all_exon_fisher <- function(){
  all_sums_list = vector("list", length(all_pairs.his.sig_list))
  for (i in 1: length(all_pairs.his.sig_list)){
    # for (i in 1:1){
    print(paste("His: ", i, sep=''))
    all_pairs.his.sig = all_pairs.his.sig_list[[i]]
    combine_exp_his = as.data.table(all_pairs.exp.padj <= 0.05 & all_pairs.his.sig < 1)
    combine_exp_his[combine_exp_his == TRUE] = 1.0
    combine_exp_his[combine_exp_his == FALSE] = 0.0 
    combine_exp_his$gene_id = all_pairs.exp.padj$gene_id
    combine_exp_his$exon_id = all_pairs.exp.padj$exon_id
    combine_exp_his = combine_exp_his %>%
      dplyr::group_split(gene_id)
    combine_exp_his = lapply(combine_exp_his, function(x) exon_cor(x))
    combine_exp_his = unlist(combine_exp_his, use.names=FALSE)
    all_sums_list[[i]] = combine_exp_his
  }
  return(all_sums_list)
}

fisher.test(df[1,], df[2,])$p
exon_cor <- function(all_exons){
  all_exons = t(as.matrix(all_exons[,3:17]))
  all_exons = p_value_calculator(cor(all_exons), ncol(all_exons))
  if ( length(table(all_exons <= 0.05 & all_exons > 0 & !is.na(all_exons))) > 1)
    return(1)
  else
    return(0)
}
all_exon_cors = get_all_exon_cor()
# saveRDS(all_exon_cors, "all_exon_cors.RDS")
temp = all_exon_cors
temp = lapply(temp, as.data.table)
temp = lapply(temp, function(x) cbind(all_genes, x))
all_exons_res = vector("list", length(temp))
for (i in 1:length(temp)){
  cor_genes = temp[[i]]
  cor_genes = cor_genes[cor_genes[[2]] == 1, 1]
  # cor_genes = lapply(cor_genes, function(x) paste(x, collapse = ', '))
  all_exons_res[[i]] = cor_genes
}
all_exons_res = Reduce(union, all_exons_res)
temp = lapply(all_exons_res, function(x) gsub("\\+", ", ", x))
temp = lapply(temp, function(x)  paste(x, collapse = ', '))[[1]]
temp = strsplit(temp, split = ', ')[[1]]
exon_res_p = temp
exon_res = readRDS("exon_res_cor.RDS")
saveRDS(exon_res, "exon_res_cor.RDS")




#========COMBINE=========



temp = all_epigenes_chi_list
temp = Reduce(union, temp)
chi_res = readRDS("tissue_spec_genes/all_epigenes_chi_list.RDS")
chi_res = Reduce(union, chi_res)

library(VennDiagram)
venn.plot <- venn.diagram(x = list(
  "Gene-level (Chi-square) 1" = chi_res,
  "Gene-level (Chi-square)" = all_chis,
  "Exon-level " = exon_res,
  "Exon-level p-value" = exon_res_p
), 
filename=NULL, col="black", fill = c("yellow", "red", "purple", "pink"), margin=0.1, alpha = 0.5,
cat.dist = 0.1, print.mode = "raw")
jpeg("overlap_gene_exon2.jpg", res = 300, width = 7, height = 5, units="in"); 
grid.draw(venn.plot); 
dev.off()










