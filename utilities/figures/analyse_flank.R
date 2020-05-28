library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(rtracklayer)
library(Hmisc)


p_value_calculator <- function(r, nrow){
  P <- r*sqrt(nrow-2)/sqrt(1-r*r)
  P <- 2*pt(-abs(P), nrow-2)
  return(P)
}
OR_func <- function(a, b){
  return(a | b)
}
get_tissue_specific_gene <- function(all_epigenes, pair_list){
  all_epigenes_intersect = Reduce(intersect, all_epigenes[pair_list])
  all_epigenes_intersect = gsub("\\+", ", ", all_epigenes_intersect)
  all_epigenes_intersect = paste(all_epigenes_intersect, collapse = ", ")
}

esc = c(1:5) #Esc
mes = c(1,6,7,8,9) #Mesendoderm
tro = c(2,6,10,11,12) #Trophoblast
gas = c(3,7,10,13,14) #Gastric
ven = c(4,8,11,13,15) #Ventricle
mus = c(5,9,12,14,15) #Muscles
tissue_type_list = list(esc, mes, tro, gas, ven, mus)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")
all_pairs.exp = list.files("exp", full.names = TRUE)
all_pairs.his = list.files("subset3", full.names = TRUE)

# exp1 = fread(all_pairs.exp[[1]])
# fwrite(exp1[,1:2], "exp_id.txt", sep = '\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

get_all_pairs.exp <- function(all_pairs.exp){
  pair.exp_list = vector("list", length(all_pairs.exp))
  for (i in 1:length(all_pairs.exp)){
    pair.exp = all_pairs.exp[[i]]
    pair.exp_list[[i]] = fread(pair.exp)$padj
  }
  pair.exp_list = lapply(pair.exp_list, function(x) ifelse((x < 0.05 & !is.na(x)), TRUE, FALSE))
  pair.exp_list = as.data.table(pair.exp_list)
  exp_id = fread("exp_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.exp_list = cbind(exp_id, pair.exp_list)
  pair.exp_list = pair.exp_list[order(pair.exp_list$V1)]
  return(pair.exp_list)
}
all_pairs.exp.padj = get_all_pairs.exp(all_pairs.exp)
colnames(all_pairs.exp.padj) = c("gene_id", "exon_id", seq(1,(ncol(all_pairs.exp.padj)-2),1))
head(all_pairs.exp.padj)
all_genes = as.data.table(unique(all_pairs.exp.padj$gene_id))
fwrite(all_genes, "gene_id.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
all_genes = fread("gene_id.txt", header = FALSE)

# flank = import.gff("reference_genome.fl300.gtf")
# flank_id = unique(as.data.table(cbind(flank$gene, flank$exon)))
# fwrite(flank_id, "flank_id.txt", sep = '\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

get_all_pairs.his <- function(all_pairs.his){
  pair.his_list = vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)){
    print(paste("Pair: ", i, sep=''))
    pair.his = all_pairs.his[[i]]
    pair.his = fread(pair.his)
    pair.his = pair.his %>%
      mutate(m_val = as.numeric(as.character(V10)),
             p_val = as.numeric(as.character(V11)),
             sig = if_else( abs(m_val) >= 1 & p_val <= 0.05 & !is.na(p_val), TRUE, FALSE)) %>%
      select(sig)
    pair.his_list[[i]] = pair.his
  }
  pair.his_list = as.data.table(pair.his_list)
  pair.his_list = pair.his_list %>%
    group_by(group = gl(n()/2, 2)) %>%
    summarise_all(function(x) Reduce(OR_func, x) ) %>%
    select(-group)
  his_id = fread("flank_id.txt", sep = '\t', quote=FALSE, header = FALSE)
  pair.his_list = cbind(his_id, pair.his_list)
  pair.his_list = pair.his_list[order(pair.his_list$V1)]
  return(pair.his_list)
}
all_pairs.his.sig = get_all_pairs.his(all_pairs.his)
colnames(all_pairs.his.sig) = c("gene_id", "exon_id", seq(1, (ncol(all_pairs.his.sig)-2),1))
head(all_pairs.his.sig)

table(all_pairs.his.sig$gene_id == all_pairs.exp.padj$gene_id)
table(all_pairs.his.sig > 0)

#============================Correlation============================
new_cor <- function(exp, his){
  if (length(unique(exp)) > 1 & length(unique(his)) > 1){
    return(cor(exp, his))
  }
  else {
    return(0.0)
  }
}
get_correlation <- function(){
  all_cors_list = vector("list", ncol(all_pairs.exp.padj)-2)
  for (i in 1:length(all_pairs.exp)){
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp.padj[[i+2]]
    his = all_pairs.his.sig[[i+2]]
    cor_table = as.data.table(cbind(exp, his))
    cor_table = cor_table %>%
      group_by(all_pairs.exp.padj$gene_id) %>%
      summarise(r_val = new_cor(exp, his),
                p_val = p_value_calculator(r_val, n()),
                p_adj = p.adjust(p_val, "fdr", n())) %>%
      select(r_val, p_val, p_adj)
    all_cors_list[[i]] = cor_table
  }
  all_cors_list = as.data.table(all_cors_list)
  all_cors_list = cbind(all_genes, all_cors_list)
  return(all_cors_list)
}
all_cors = get_correlation()

all_p_val = as.data.frame(all_cors)[, c(1, seq(0, ncol(all_cors), 3))]
all_p_adj = as.data.frame(all_cors)[, seq(1, ncol(all_cors), 3)]
all_r_val = as.data.frame(all_cors)[, c(1, seq(2, ncol(all_cors), 3))]

all_epigenes_cor = vector("list", ncol(all_r_val)-1)
for (i in 2:ncol(all_p_adj)){
  all_epigenes_cor[[i-1]] = all_p_adj[all_p_adj[[i]] <= 0.05  & !is.na(all_p_adj[[i]]), 1] 
}
saveRDS(all_epigenes_cor, "all_epigenes_cor.RDS")

all_epigenes_intersect = Reduce(intersect, all_epigenes_cor)
all_epigenes_intersect = gsub("\\+", ", ", all_epigenes_intersect)
paste(all_epigenes_intersect, collapse = ", ")

all_epigenes_union = Reduce(union, all_epigenes_cor)
all_epigenes_union = gsub("\\+", ", ", all_epigenes_union)
paste(all_epigenes_union, collapse = ", ")

lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_cor, x))

#============================Fishers Exact============================
new_fisher.test <- function(exp, his){
  exp = as.factor(exp)
  his = as.factor(his)
  if ( nlevels(exp) > 1 & nlevels(his) > 1){
    return( fisher.test(exp, his)$p.value )
  }
  else {
    return(1.0)
  }
}
get_fisher <- function(){
  all_fishers_list = vector("list", ncol(all_pairs.exp.padj)-2)
  for (i in 1:length(all_pairs.exp)){
  # for (i in 1:1){
    print(paste("Pair: ", i, sep=''))
    exp = all_pairs.exp.padj[[i+2]]
    his = all_pairs.his.sig[[i+2]]
    fisher_table = as.data.table(cbind(exp, his))
    fisher_table = fisher_table %>%
      group_by(all_pairs.exp.padj$gene_id) %>%
      summarise(p_val = new_fisher.test(exp, his)) %>%
      select(p_val)
    all_fishers_list[[i]] = fisher_table
  }
  all_fishers_list = as.data.table(all_fishers_list)
  all_fishers_list = cbind(all_genes, all_fishers_list)
  return(all_fishers_list)
}
all_fishers = get_fisher()

all_epigenes_fisher = vector("list", ncol(all_fishers)-1)
for (i in 2:ncol(all_fishers)){
  all_epigenes_fisher[i-1] = all_fishers[all_fishers[[i]] <= 0.05  & !is.na(all_fishers[[i]]), 1]
}
saveRDS(all_epigenes_fisher, "all_epigenes_fisher.RDS")

all_epigenes_fisher_intersect = Reduce(intersect, all_epigenes_fisher)
all_epigenes_fisher_intersect = gsub("\\+", ", ", all_epigenes_fisher_intersect)
paste(all_epigenes_fisher_intersect, collapse = ", ")

all_epigenes_fisher_union = Reduce(union, all_epigenes_fisher)
all_epigenes_fisher_union = gsub("\\+", "\\, ", all_epigenes_fisher_union)
paste(all_epigenes_fisher_union, collapse = ", ")

lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_fisher, x))

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
get_chi <- function(){
  all_chis_list = vector("list", ncol(all_pairs.exp.padj)-2)
  for (i in 1:length(all_pairs.exp)){
    # for (i in 1:1){
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
all_chis = get_chi()

all_epigenes_chi = vector("list", ncol(all_chis)-1)
for (i in 2:ncol(all_chis)){
  all_epigenes_chi[i-1] = all_chis[all_chis[[i]] <= 0.05  & !is.na(all_chis[[i]]), 1]
}
saveRDS(all_epigenes_chi, "all_epigenes_chi.RDS")

all_epigenes_chi_intersect = Reduce(intersect, all_epigenes_chi)
all_epigenes_chi_intersect = gsub("\\+", ", ", all_epigenes_chi_intersect)
paste(all_epigenes_chi_intersect, collapse = ", ")

all_epigenes_chi_union = Reduce(union, all_epigenes_chi)
all_epigenes_chi_union = gsub("\\+", "\\, ", all_epigenes_chi_union)
paste(all_epigenes_chi_union, collapse = ", ")

lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_epigenes_chi, x))


#===========================COMBINE===========================
paste(Reduce(intersect, list(all_epigenes_union, all_epigenes_fisher_union, all_epigenes_chi_union)), collapse=', ')

library("VennDiagram", quietly=TRUE)
venn.plot <- venn.diagram(x = list(
  "From Correlation" = all_epigenes_union,
  "From Chi-square" = all_epigenes_chi_union,
  "From Fisher's Exact Test" = all_epigenes_fisher_union
), 
filename=NULL, col="black", fill = c("yellow", "red", "purple"), margin=0.1, alpha = 0.5,
cat.dist = 0.1, lwd = 0.5)
jpeg("compare_test_fl100.jpg", res = 300, width = 7, height = 5, units="in"); 
grid.draw(venn.plot); 
dev.off()

