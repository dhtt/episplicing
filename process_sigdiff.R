library(data.table, quietly=TRUE)
library(stringr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(org.Hs.eg.db)
library(GOSim)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine")
gene_id = data.frame(fread("gene_id", header=FALSE, col.names=c("gene_id")))
epi_id1 = "E003"
epi_id2 = "E004"
epi_id3 = "E005"
prepare_metdiff <- function(metdiff_tables_list){
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table$V14 = as.numeric(as.character(metdiff_table$V14))
    metdiff_table$V15 = as.numeric(as.character(metdiff_table$V15))
    metdiff_table$V16 = as.numeric(as.character(metdiff_table$V16))
    metdiff_table = metdiff_table[, c(1,4,5,7,9,14,15,16)]
    colnames(metdiff_table) = c("chr", "start", "end", "strand", "id", "p-val", "q-val", "metdiff")
    metdiff_tables_list[[i]] = metdiff_table[metdiff_table$`q-val` <= 0.05 & !is.na(metdiff_table$`q-val`) & abs(metdiff_table$`metdiff`) >= 75, ]
  }
  return(metdiff_tables_list)
}
get_sig_genes_meth <- function(metdiff_tables_list, mode){
  sig_genes_meth_list = vector("list", 2)
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    gene_id = as.data.table(str_split_fixed(metdiff_table$`id`, "\\;|\\ ", 8))[,2]
    gene_id <- gene_id %>%
      mutate(gene_id = gsub("\"", "", V2)) %>%
      dplyr::select(gene_id) %>%
      unique()
    sig_genes_meth_list[[i]] = c(gene_id[,1])
  }
  if (mode == "none"){
    return(sig_genes_meth_list)
  } 
  else if (mode == "intersect") {
    return(intersect(sig_genes_meth_list[[1]],sig_genes_meth_list[[2]]))
  }
  else if (mode == "union") {
    return(union(sig_genes_meth_list[[1]],sig_genes_meth_list[[2]]))
  }
}
prepare_hisdiff <- function(hisdiff_tables_list){
  for (i in range(1:2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table$V13 = as.numeric(as.character(hisdiff_table$V13))
    hisdiff_table$V14 = as.numeric(as.character(hisdiff_table$V14))
    hisdiff_table = hisdiff_table[, c(1,4,5,7,9,13,14)]
    colnames(hisdiff_table) = c("chr", "start", "end", "strand", "id", "m-val", "p-val")
    hisdiff_tables_list[[i]] = hisdiff_table[!is.na(hisdiff_table$`p-val`) & !is.na(hisdiff_table$`m-val`) & 
                                             hisdiff_table$`p-val` <= 0.05 & abs(hisdiff_table$`m-val`) >= 1,]
  }
  return(hisdiff_tables_list)
}
get_sig_genes_his <- function(hisdiff_tables_list, mode){
  sig_genes_his_list = vector("list", 2)
  for (i in range(1:2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    gene_id = as.data.table(str_split_fixed(hisdiff_table$`id`, "\\;|\\ ", 8))[,2]
    gene_id <- gene_id %>%
      mutate(gene_id = gsub("\"", "", V2)) %>%
      dplyr::select(gene_id) %>%
      unique()
    sig_genes_his_list[[i]] = c(gene_id[,1])
  }
  if (mode == "none"){
    return(sig_genes_his_list)
  } 
  else if (mode == "intersect") {
    return(intersect(sig_genes_his_list[[1]], sig_genes_his_list[[2]]))
  }
  else if (mode == "union") {
    return(union(sig_genes_his_list[[1]], sig_genes_his_list[[2]]))
  }
}
prepare_expdiff <- function(expdiff_tables_list){
  for (i in range(1:length(expdiff_tables_list))){
    expdiff_table = expdiff_tables_list[[i]]
    expdiff_table <- as.data.frame(expdiff_table[expdiff_table$`padj` <= 0.05 & !is.na(expdiff_table$`padj`),])
    expdiff_tables_list[[i]] = expdiff_table
    }
  return(expdiff_tables_list)
}
get_sig_genes_exp <- function(expdiff_tables_list, mode){
  #TODO also intron
  sig_genes_exp_list = vector("list", length(expdiff_tables_list))
  for (i in range(1:length(expdiff_tables_list))){
    expdiff_table = expdiff_tables_list[[i]]
    gene_id <- expdiff_table %>%
      dplyr::select(groupID) %>%
      unique()
    sig_genes_exp_list[[i]] = c(gene_id[,1])
  }
  if (mode == "none"){
    return(sig_genes_exp_list[[1]])
  } 
  # else if (mode == "intersect") {
  #   return(intersect(sig_genes_exp_list[[1]], sig_genes_exp_list[[2]]))
  # }
  # else if (mode == "union") {
  #   return(union(sig_genes_exp_list[[1]], sig_genes_exp_list[[2]]))
  # }
}

#===== Pair 1 =====
#===== Meth Diff =====
file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
pair = paste(epi_id1, epi_id2, sep='_')
met.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
met.diff.exon1 = fread(met.diff1[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
met.diff.pi1 = fread(met.diff1[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
metdiff_tables_list1 = list(met.diff.exon1, met.diff.pi1)
metdiff_tables_list1 = prepare_metdiff(metdiff_tables_list1)
sig_genes_meth_list1 = get_sig_genes_meth(metdiff_tables_list1, mode="union")

#===== H3K36me3 Diff =====
file_path = normalizePath(paste(getwd(), "H3K36me3", "annotatedcounts", sep='/'))
pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
             paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
H3K36me3.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
H3K36me3.exon1 = fread(H3K36me3.diff1[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3.pi1 = fread(H3K36me3.diff1[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3diff_tables_list1 = list(H3K36me3.exon1, H3K36me3.pi1)
H3K36me3diff_tables_list1 = prepare_hisdiff(H3K36me3diff_tables_list1)
sig_genes_his_list1 = get_sig_genes_his(H3K36me3diff_tables_list1, mode="union")

#===== Expression Diff =====
file_path = normalizePath(paste(getwd(), "expression", sep='/'))
pair = paste(epi_id1, epi_id2, "res.csv", sep='_')
exp.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
exp.diff.exon1 = fread(exp.diff1[[1]], header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
expdiff_tables_list1 = list(exp.diff.exon1)
expdiff_tables_list1 = prepare_expdiff(expdiff_tables_list1)
sig_genes_exp_list1 = get_sig_genes_exp(expdiff_tables_list1, mode = "none")

#===== Differential genes =====
E003_E004_expmet = intersect(sig_genes_exp_list1, sig_genes_meth_list1)
E003_E004_exphis = intersect(sig_genes_exp_list1, sig_genes_his_list1)
E003_E004_all = intersect(intersect(sig_genes_exp_list1, sig_genes_meth_list1), sig_genes_his_list1)
E003_E004_union = union(E003_E004_expmet, E003_E004_exphis)

#===== Pair 2 =====
#===== Meth Diff =====
file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
pair = paste(epi_id1, epi_id3, sep='_')
met.diff2 = list.files(file_path, pattern=pair, full.names=TRUE)
met.diff.exon2 = fread(met.diff2[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
met.diff.pi2 = fread(met.diff2[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
metdiff_tables_list2 = list(met.diff.exon2, met.diff.pi2)
metdiff_tables_list2 = prepare_metdiff(metdiff_tables_list2)
sig_genes_meth_list2 = get_sig_genes_meth(metdiff_tables_list2, mode="union")

#===== H3K36me3 Diff =====
file_path = normalizePath(paste(getwd(), "H3K36me3", "annotatedcounts", sep='/'))
pair = paste(paste('^pi', paste(epi_id1, epi_id3, sep='_'), sep='_'),
             paste('^exon', paste(epi_id1, epi_id3, sep='_'), sep='_'), sep='|')
H3K36me3.diff2 = list.files(file_path, pattern=pair, full.names=TRUE)
H3K36me3.exon2 = fread(H3K36me3.diff2[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3.pi2 = fread(H3K36me3.diff2[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3diff_tables_list2 = list(H3K36me3.exon2, H3K36me3.pi2)
H3K36me3diff_tables_list2 = prepare_hisdiff(H3K36me3diff_tables_list2)
sig_genes_his_list2 = get_sig_genes_his(H3K36me3diff_tables_list2, mode="union")

#===== Expression Diff =====
file_path = normalizePath(paste(getwd(), "expression", sep='/'))
pair = paste(epi_id1, epi_id3, "res.csv", sep='_')
exp.diff2 = list.files(file_path, pattern=pair, full.names=TRUE)
exp.diff.exon2 = fread(exp.diff2[[1]], header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
expdiff_tables_list2 = list(exp.diff.exon2)
expdiff_tables_list2 = prepare_expdiff(expdiff_tables_list2)
sig_genes_exp_list2 = get_sig_genes_exp(expdiff_tables_list2, mode = "none")

#===== Differential genes =====
E003_E005_expmet = intersect(sig_genes_exp_list2, sig_genes_meth_list2)
E003_E005_exphis = intersect(sig_genes_exp_list2, sig_genes_his_list2)
E003_E005_all = intersect(intersect(sig_genes_exp_list2, sig_genes_meth_list2), sig_genes_his_list2)
E003_E005_union = union(E003_E005_expmet, E003_E005_exphis)

#===== Pair 3 =====
#===== Meth Diff =====
file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
pair = paste(epi_id2, epi_id3, sep='_')
met.diff3 = list.files(file_path, pattern=pair, full.names=TRUE)
met.diff.exon3 = fread(met.diff3[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
met.diff.pi3 = fread(met.diff3[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
metdiff_tables_list3 = list(met.diff.exon3, met.diff.pi3)
metdiff_tables_list3 = prepare_metdiff(metdiff_tables_list3)
sig_genes_meth_list3 = get_sig_genes_meth(metdiff_tables_list3, mode="union")

#===== H3K36me3 Diff =====
file_path = normalizePath(paste(getwd(), "H3K36me3", "annotatedcounts", sep='/'))
pair = paste(paste('^pi', paste(epi_id2, epi_id3, sep='_'), sep='_'),
             paste('^exon', paste(epi_id2, epi_id3, sep='_'), sep='_'), sep='|')
H3K36me3.diff3 = list.files(file_path, pattern=pair, full.names=TRUE)
H3K36me3.exon3 = fread(H3K36me3.diff3[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3.pi3 = fread(H3K36me3.diff3[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
H3K36me3diff_tables_list3 = list(H3K36me3.exon3, H3K36me3.pi3)
H3K36me3diff_tables_list3 = prepare_hisdiff(H3K36me3diff_tables_list3)
sig_genes_his_list3 = get_sig_genes_his(H3K36me3diff_tables_list3, mode="union")

#===== Expression Diff =====
file_path = normalizePath(paste(getwd(), "expression", sep='/'))
pair = paste(epi_id2, epi_id3, "res.csv", sep='_')
exp.diff3 = list.files(file_path, pattern=pair, full.names=TRUE)
exp.diff.exon3 = fread(exp.diff3[[1]], header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
expdiff_tables_list3 = list(exp.diff.exon3)
expdiff_tables_list3 = prepare_expdiff(expdiff_tables_list3)
sig_genes_exp_list3 = get_sig_genes_exp(expdiff_tables_list3, mode = "none")

#===== Differential genes =====
E004_E005_expmet = intersect(sig_genes_exp_list3, sig_genes_meth_list3)
E004_E005_exphis = intersect(sig_genes_exp_list3, sig_genes_his_list3)
E004_E005_all = intersect(intersect(sig_genes_exp_list3, sig_genes_meth_list3), sig_genes_his_list3)
E004_E005_union = union(E004_E005_expmet, E004_E005_exphis)

#===== ANNOTATION =====
library(topGO)
library(goProfiles)
all_genes = unique(exp.diff.exon1$groupID)

get_entrez_id <-function(gene_list){
  entrez_gene_list = select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}
E003_E004_expmet.entrez = get_entrez_id(E003_E004_expmet)
E003_E004_exphis.entrez = get_entrez_id(E003_E004_exphis)
E003_E004_all.entrez = get_entrez_id(E003_E004_all)

E003_E005_expmet.entrez = get_entrez_id(E003_E005_expmet)
E003_E005_exphis.entrez = get_entrez_id(E003_E005_exphis)
E003_E005_all.entrez = get_entrez_id(E003_E005_all)

E004_E005_expmet.entrez = get_entrez_id(E004_E005_expmet)
E004_E005_exphis.entrez = get_entrez_id(E004_E005_exphis)
E004_E005_all.entrez = get_entrez_id(E004_E005_all)
# 
# E003_E004_all_GO = getGOInfo(E003_E004_all.entrez)
# E003_E005_all_GO = getGOInfo(E003_E005_all.entrez)
# E004_E005_all_GO = getGOInfo(E004_E005_all.entrez)

all_pairs = intersect(E004_E005_all, intersect(E003_E004_all, E003_E005_all))
all_pairs.entrez = get_entrez_id(all_pairs)
all_genes.entrez = get_entrez_id(all_genes)

paste(E003_E004_all, collapse = ',')



temp1 = basicProfile(E003_E004_all.entrez[!is.na(E003_E004_all.entrez)], idType = "Entrez", onto="BP", level = 2, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = TRUE)
temp2 = basicProfile(E003_E005_all.entrez[!is.na(E003_E005_all.entrez)], idType = "Entrez", onto="BP", level = 2, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = TRUE)
temp3 = basicProfile(E004_E005_all.entrez[!is.na(E004_E005_all.entrez)], idType = "Entrez", onto="BP", level = 2, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = TRUE)
temp4 = basicProfile(all_pairs.entrez, idType = "Entrez", onto="MF", level = 2, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = TRUE)
merge_temp = mergeProfilesLists(temp1, temp2, profNames = c("H1ESC_Mesendoderm", "H1ESC_Trophoblast"))
merge_temp = mergeProfilesLists(temp1, temp3, profNames = c("H1ESC_Mesendoderm", "H1ESC_Trophoblast"))

plotProfiles(c(temp1, temp2, temp3), percentage=TRUE, multiplePlots = FALSE, HORIZVERT = TRUE,
             legendText = TRUE)

temp1$MF



