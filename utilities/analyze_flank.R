library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(rtracklayer)
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
# all_genes = as.data.table(unique(all_pairs.exp.padj$gene_id))
# fwrite(all_genes, "gene_id.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
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
# all_pairs.his = list.files(histone_type_list[[3]], full.names = TRUE)
# all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
# colnames(all_pairs.his.sig) = c("gene_id", "exon_id", seq(1, (ncol(all_pairs.his.sig)-2),1))
# head(all_pairs.his.sig)
temp = all_pairs.his.sig_list[[6]]
table(all_pairs.his.sig$gene_id == all_pairs.exp.padj$gene_id)

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

temp = join_all_his_res(all_epigenes_cor_pearson_list)
temp = lapply(tissue_type_list, function(x) get_tissue_specific_gene(temp, x))
lapply(temp, function(x) length(strsplit(x[1], split = ', ')[[1]]))
temp = transpose(as.data.table(temp))
fwrite(temp, "cor_pearson_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')

################################################
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

temp = join_all_his_res(all_epigenes_cor_spearman_list)
temp = lapply(tissue_type_list, function(x) get_tissue_specific_gene(temp, x))
lapply(temp, function(x) length(strsplit(x[1], split = ', ')[[1]]))
temp = transpose(as.data.table(temp))
fwrite(temp, "cor_spearman_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')
##============================Get Chi============================
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
      select(p_val)
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

temp = join_all_his_res(all_epigenes_chi_list)
temp = lapply(tissue_type_list, function(x) get_tissue_specific_gene(temp, x))
lapply(temp, function(x) length(strsplit(x[1], split = ', ')[[1]]))
temp = transpose(as.data.table(temp))
fwrite(temp, "chi_tissue_spec.txt", quote = FALSE, row.names = TRUE, sep='\t')

#===========================COMBINE===========================
paste(Reduce(intersect, list(all_epigenes_union_pearson, all_epigenes_union_spearman, all_epigenes_chi_union)), collapse=', ')

library("VennDiagram", quietly=TRUE)
venn.plot <- venn.diagram(x = list(
  "From Pearsons Correlation" = all_epigenes_union_pearson,
  "From Spearmans Correlation" = all_epigenes_union_spearman,
  "From Chi-square" = all_epigenes_chi_union
), 
filename=NULL, col="black", fill = c("yellow", "red", "purple"), margin=0.1, alpha = 0.5,
cat.dist = 0.1, lwd = 0.5)
jpeg("compare_test_fl200_cat.jpg", res = 300, width = 7, height = 5, units="in"); 
grid.draw(venn.plot); 
dev.off()

