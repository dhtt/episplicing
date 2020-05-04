library(data.table, quietly=TRUE)
library(stringr, quietly=TRUE)
library(dplyr, quietly=TRUE)

setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine2")
gene_id = data.frame(fread("gene_id", header=FALSE, col.names=c("gene_id")))
epi_id1 = "E003"
epi_id2 = "E004"
epi_id3 = "E005"
epi_id4 = "E094"
epi_id5 = "E095"
epi_id6 = "E096"
prepare_metdiff <- function(metdiff_tables_list){
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table$V14 = as.numeric(as.character(metdiff_table$V14))
    metdiff_table$V15 = as.numeric(as.character(metdiff_table$V15))
    metdiff_table$V16 = as.numeric(as.character(metdiff_table$V16))
    metdiff_table = metdiff_table[, c(1,4,5,7,9,14,15,16)]
    colnames(metdiff_table) = c("chr", "start", "end", "strand", "id", "p-val", "q-val", "metdiff")
    metdiff_tables_list[[i]] = metdiff_table[metdiff_table$`q-val` <= 0.05 & !is.na(metdiff_table$`q-val`) & abs(metdiff_table$`metdiff`) >= 25, ]
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
get_met_diff <- function(epi_id1, epi_id2, mode="union"){
  file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  met.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
  met.diff.exon1 = fread(met.diff1[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  met.diff.pi1 = fread(met.diff1[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_tables_list1 = list(met.diff.exon1, met.diff.pi1)
  metdiff_tables_list1 = prepare_metdiff(metdiff_tables_list1)
  sig_genes_meth_list1 = get_sig_genes_meth(metdiff_tables_list1, mode=mode)
  return(sig_genes_meth_list1)
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
get_his_diff <- function(type, epi_id1, epi_id2, mode="union"){
  file_path = normalizePath(paste(getwd(), type, "annotatedcounts", sep='/'))
  pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
               paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
  H3K36me3.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
  H3K36me3.exon1 = fread(H3K36me3.diff1[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  H3K36me3.pi1 = fread(H3K36me3.diff1[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  H3K36me3diff_tables_list1 = list(H3K36me3.exon1, H3K36me3.pi1)
  H3K36me3diff_tables_list1 = prepare_hisdiff(H3K36me3diff_tables_list1)
  sig_genes_his_list1 = get_sig_genes_his(H3K36me3diff_tables_list1, mode=mode)
  return(sig_genes_his_list1)
}

get_exp_diff <- function(epi_id1, epi_id2, mode="none"){
  file_path = normalizePath(paste(getwd(), "expression","annotatedcounts", sep='/'))
  pair = paste(epi_id1, epi_id2, "res.csv", sep='_')
  exp.diff1 = list.files(file_path, pattern=pair, full.names=TRUE)
  exp.diff.exon1 = fread(exp.diff1[[1]], header=TRUE, stringsAsFactors = FALSE, quote=FALSE)
  expdiff_tables_list1 = list(exp.diff.exon1)
  expdiff_tables_list1 = prepare_expdiff(expdiff_tables_list1)
  sig_genes_exp_list1 = get_sig_genes_exp(expdiff_tables_list1, mode = mode)
  return(sig_genes_exp_list1)
}
prepare_expdiff <- function(expdiff_tables_list){
  for (i in range(1:length(expdiff_tables_list))){
    expdiff_table = expdiff_tables_list[[i]]
    expdiff_table <- as.data.frame(expdiff_table[expdiff_table$`padj` <= 0.05 & !is.na(expdiff_table$`padj`),])
    expdiff_tables_list[[i]] = expdiff_table
    }
  return(expdiff_tables_list)
}
get_sig_genes_exp <- function(expdiff_tables_list, mode="union"){
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
pair1_2.met = get_met_diff(epi_id1, epi_id2)
pair1_2.h3k36me3 = get_his_diff("H3K36me3",epi_id1, epi_id2) 
pair1_2.h3k27ac = get_his_diff("H3K27ac",epi_id1, epi_id2) 
pair1_2.exp = get_exp_diff(epi_id1, epi_id2) 

pair1_2.epg = intersect(pair1_2.exp, union(pair1_2.met, union(pair1_2.h3k36me3, pair1_2.h3k27ac)))
paste(pair1_2.epg, collapse = ",")
#===== Pair 2 =====
pair1_3.met = get_met_diff(epi_id1, epi_id3)
pair1_3.h3k36me3 = get_his_diff("H3K36me3",epi_id1, epi_id3) 
pair1_3.h3k27ac = get_his_diff("H3K27ac",epi_id1, epi_id3) 
pair1_3.exp = get_exp_diff(epi_id1, epi_id3) 

pair1_3.epg = intersect(pair1_3.exp, union(pair1_3.met, union(pair1_3.h3k36me3, pair1_3.h3k27ac)))
paste(pair1_3.epg, collapse = ",")
#===== Pair 3 =====
pair2_3.met = get_met_diff(epi_id2, epi_id3)
pair2_3.h3k36me3 = get_his_diff("H3K36me3",epi_id2, epi_id3) 
pair2_3.h3k27ac = get_his_diff("H3K27ac",epi_id2, epi_id3) 
pair2_3.exp = get_exp_diff(epi_id2, epi_id3) 

pair2_3.epg = intersect(pair2_3.exp, union(pair2_3.met, union(pair2_3.h3k36me3, pair2_3.h3k27ac)))
paste(pair2_3.epg, collapse = ",")
#===== Pair 4 =====
pair4_5.met = get_met_diff(epi_id4, epi_id5)
pair4_5.h3k36me3 = get_his_diff("H3K36me3",epi_id4, epi_id5) 
pair4_5.h3k27ac = get_his_diff("H3K27ac",epi_id4, epi_id5) 
pair4_5.exp = get_exp_diff(epi_id4, epi_id5) 

pair4_5.epg = intersect(pair4_5.exp, union(pair4_5.met, union(pair4_5.h3k36me3, pair4_5.h3k27ac)))
paste(pair4_5.epg, collapse = ",")

#===== Pair 5 =====
pair4_6.met = get_met_diff(epi_id4, epi_id6)
pair4_6.h3k36me3 = get_his_diff("H3K36me3", epi_id4, epi_id6) 
pair4_6.h3k27ac = get_his_diff("H3K27ac", epi_id4, epi_id6) 
pair4_6.exp = get_exp_diff(epi_id4, epi_id6) 

pair4_6.epg = intersect(pair4_6.exp, union(pair4_6.met, union(pair4_6.h3k36me3, pair4_6.h3k27ac)))
paste(pair2_3.epg, collapse = ",")

#===== Pair 6 =====
pair5_6.met = get_met_diff(epi_id5, epi_id6)
pair5_6.h3k36me3 = get_his_diff("H3K36me3", epi_id5, epi_id6) 
pair5_6.h3k27ac = get_his_diff("H3K27ac", epi_id5, epi_id6) 
pair5_6.exp = get_exp_diff(epi_id5, epi_id6) 

pair5_6.epg = intersect(pair5_6.exp, union(pair5_6.met, union(pair5_6.h3k36me3, pair5_6.h3k27ac)))
paste(pair2_3.epg, collapse = ",")

#===== Differential genes =====
tot_epg = list(pair1_2.epg, pair1_3.epg, pair2_3.epg, pair4_5.epg, pair4_6.epg, pair5_6.epg)
tot_epg = as.data.frame(Reduce(intersect, tot_epg))
grep("NANOG", pair1_2.exp)
grep("NANOG", pair1_2.h3k27ac)
grep("NANOG", pair1_2.h3k36me3)
grep("NANOG", pair1_2.met)

grep("SRP19", pair4_5.exp)
grep("SRP19", pair4_5.h3k27ac)
grep("SRP19", pair4_5.h3k36me3)
grep("SRP19", pair4_5.met)

#===== Summary  =====
tot_met = list(pair1_2.met, pair1_3.met, pair2_3.met, pair4_5.met, pair4_6.met, pair5_6.met)
tot_h3k36me3 = list(pair1_2.h3k36me3, pair1_3.h3k36me3, pair2_3.h3k36me3,
                    pair4_5.h3k36me3, pair4_6.h3k36me3, pair5_6.h3k36me3)
tot_h3k27ac = list(pair1_2.h3k27ac, pair1_3.h3k27ac, pair2_3.h3k27ac,
                    pair4_5.h3k27ac, pair4_6.h3k27ac, pair5_6.h3k27ac)
tot_exp = list(pair1_2.exp, pair1_3.exp, pair2_3.exp, pair4_5.exp, pair4_6.exp, pair5_6.exp)
#-------
all_met_genes = Reduce(union, tot_met)
all_h3k36me3_genes = Reduce(union, tot_h3k36me3)
all_h3k27ac_genes = Reduce(union, tot_h3k27ac)
all_exp_genes = Reduce(union, tot_exp)
# 
# example_genes = as.data.frame(intersect(all_exp_genes, intersect(all_h3k27ac_genes, intersect(all_met_genes, all_h3k36me3_genes))))
# head(example_genes)
# ex1 = "SMAD3"
# ex1 = "NANOG;"

#-------
exp_by_tissue = list(union(pair1_2.exp, pair1_3.exp), union(pair1_2.exp, pair2_3.exp),
     union(pair1_3.exp, pair2_3.exp), union(pair4_5.exp, pair4_6.exp),
     union(pair4_5.exp, pair5_6.exp), union(pair4_6.exp, pair5_6.exp))
exp_by_tissue.no = as.data.table(lapply(exp_by_tissue, function(x) length(x)))

met_by_tissue = list(union(pair1_2.met, pair1_3.met), union(pair1_2.met, pair2_3.met),
                     union(pair1_3.met, pair2_3.met), union(pair4_5.met, pair4_6.met),
                     union(pair4_5.met, pair5_6.met), union(pair4_6.met, pair5_6.met))
met_by_tissue.no = as.data.table(lapply(met_by_tissue, function(x) length(x)))

h3k36me3_by_tissue = list(union(pair1_2.h3k36me3, pair1_3.h3k36me3), union(pair1_2.h3k36me3, pair2_3.h3k36me3),
                     union(pair1_3.h3k36me3, pair2_3.h3k36me3), union(pair4_5.h3k36me3, pair4_6.h3k36me3),
                     union(pair4_5.h3k36me3, pair5_6.h3k36me3), union(pair4_6.h3k36me3, pair5_6.h3k36me3))
h3k36me3_by_tissue.no = as.data.table(lapply(h3k36me3_by_tissue, function(x) length(x)))

h3k27ac_by_tissue = list(union(pair1_2.h3k27ac, pair1_3.h3k27ac), union(pair1_2.h3k27ac, pair2_3.h3k27ac),
                          union(pair1_3.h3k27ac, pair2_3.h3k27ac), union(pair4_5.h3k27ac, pair4_6.h3k27ac),
                          union(pair4_5.h3k27ac, pair5_6.h3k27ac), union(pair4_6.h3k27ac, pair5_6.h3k27ac))
h3k27ac_by_tissue.no = as.data.table(lapply(h3k27ac_by_tissue, function(x) length(x)))

#-------
tissue_stat = as.data.table(t(rbind(exp_by_tissue.no, met_by_tissue.no, h3k36me3_by_tissue.no, h3k27ac_by_tissue.no)))
all_tissue_stat = as.data.table(lapply(list(all_exp_genes, all_met_genes, all_h3k36me3_genes, all_h3k27ac_genes), function(x) length(x)))
summary_freq = as.data.frame(c("H1ESC", "Mesendoderm", "Trophoblast", "Gastric", "Ventricle", "Lung", "Total"))
summary_freq = cbind(summary_freq, rbind(tissue_stat, all_tissue_stat))
colnames(summary_freq) = c("Tissue","DEU", "Methylation", "H3K36me3", "H3K27ac")

#===== Table=====
library(formattable)
tissue_formatter <- formatter("span", 
                              style = x ~ style(
                                width = suffix(x, "px"),
                                font.weight = "bold", 
                                color = ifelse(x == "Total", "black", "gray")))
formattable(summary_freq,
            align =c("l","c","c","c","c"), 
            list(`Tissue` = tissue_formatter,
                 `DEU` = color_tile("white", "wheat"),
                 `Methylation` = color_tile("white", "wheat"),
                 `H3K36me3` = color_tile("white", "wheat"),
                 `H3K27ac` = color_tile("white", "wheat")
                 )
            )

#=====Venn=====
library("VennDiagram", quietly=TRUE)

venn.plot <- venn.diagram(x = list(
  "Differentially methylated genes" = all_met_genes,
  "Differentially modified H3K36me3" = all_h3k36me3_genes,
  "Differential exon usage genes" = all_exp_genes
), 
filename=NULL, col="black", fill = c("yellow", "red", "purple"), margin=0.1, alpha = 0.5)
jpeg("summary.jpg"); 
grid.draw(venn.plot); 
dev.off()

#=====Other=====

library(rcartocolor)
# l1 = E003_E004_union
# l2 = E003_E005_union
# l3 = E004_E005_union
l4 = E003_E004_union
l5 = E003_E005_union
l6 = E004_E005_union
l0 = list()
ls = list(l1,l2,l3,l4,l5,l6)
cell1 = rep(c("ESC", "Mesendoderm", "Trophoblast", "Gastric", "Ventricle", "Lung"),6)
cell2 = cell1[order(cell1)]
epi_genes = data.frame(cbind(cell2, cell1))
epi_genes$pair = paste(epi_genes$cell2, epi_genes$cell1, sep = '_')
epi_genes$no_genes = NA
all_pairs = data.frame(cbind(c("ESC", "ESC", "Mesendoderm", "Gastric", "Gastric", "Ventricle"),
                  c("Mesendoderm", "Trophoblast", "Trophoblast", "Ventricle", "Lung", "Lung")))
colnames(all_pairs) = c("cell1","cell2")
all_pairs$pair1 = paste(all_pairs$cell1, all_pairs$cell2, sep='_')
all_pairs$pair2 = paste(all_pairs$cell2, all_pairs$cell1, sep='_')
all_pairs$no_genes = c(lapply(ls, function(x) length(x)))
for (i in 1:36){
  # print(epi_genes[i,]$pair)
  for (j in 1:6){
    if (epi_genes[i,]$pair == all_pairs[j,]$pair1){
      print(all_pairs[j,]$no_genes)
      epi_genes[i,]$no_genes = all_pairs[j,]$no_genes
    }
    if (epi_genes[i,]$cell2 == epi_genes[i,]$cell1){
      epi_genes[i,]$no_genes = 0
    }
  }
}

set1 = epi_genes[grep("Gastric|Ventricle|Lung", epi_genes$cell2),]
set1[5,]$no_genes = set1[9,]$no_genes
set1[9,]$no_genes= NA
set1 = set1[grep("Gastric|Ventricle|Lung", set1$cell1),]
set2 = epi_genes[grep("ESC|Mesendoderm|Trophoblast", epi_genes$cell2),]
set2 = set2[grep("ESC|Mesendoderm|Trophoblast", set2$cell1),]
set1$no_genes = as.numeric(set1$no_genes)
set2$no_genes = as.numeric(set2$no_genes)
set1$type = "mature"
set2$type = "early_staged"
epi_genes = data_frame(rbind(set1,set2))
plot1 = ggplot(data = set1, aes(x = cell1, y = cell2)) +
  geom_tile(aes(fill = no_genes), color='white') +
  labs(x=NULL, y=NULL) +
  scale_fill_carto_c(palette = "ag_GrnYl", name = "No. of genes") +
  annotate("text", x = 1:3, y = 1:3, label = c("Gastric", "Ventricle", "Lung"), colour = "white", size=3) +
  theme(axis.text.x=element_text(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        axis.text.y=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'))
ggsave(plot = plot1, filename = "set2.png", width = 4, height = 3)
plot2 = ggplot(data = set2, aes(x = cell1, y = cell2)) +
  geom_tile(aes(fill = no_genes), color='white') +
  labs(x=NULL, y=NULL) +
  scale_fill_carto_c(palette = "ag_GrnYl", name = "No. of genes") +
  annotate("text", x = 1:3, y = 1:3, label = c("ESC", "Mesendoderm", "Trophoblast"), colour = "white", size = 3) +
  theme(axis.text.x=element_text(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'))
plot2
ggsave(plot = plot2, filename = "set1.png", width = 4, height = 3)



#===== ANNOTATION =====

intersect(E004_E005_all, intersect(E003_E004_all, E003_E005_all))

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




fwrite(list(E004_E005_union), file="004_005.u.txt", sep="\n")

temp = GOSim::GOenrichment(E003_E004_all.entrez[!is.na(E003_E004_all.entrez)], all_genes.entrez)
temp_term = temp$GOTerms
paste(temp_term$go_id, collapse = ' ')


temp1 = basicProfile(E003_E004_all.entrez[!is.na(E003_E004_all.entrez)], idType = "Entrez", onto="MF", level = 3, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = FALSE, empty.cats = FALSE)
temp2 = basicProfile(E003_E005_all.entrez[!is.na(E003_E005_all.entrez)], idType = "Entrez", onto="MF", level = 3, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = FALSE)
temp3 = basicProfile(E004_E005_all.entrez[!is.na(E004_E005_all.entrez)], idType = "Entrez", onto="MF", level = 3, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = FALSE)
temp4 = basicProfile(all_pairs.entrez, idType = "Entrez", onto="MF", level = 2, 
                     orgPackage = "org.Hs.eg.db", 
                     ord = TRUE, cat.names = TRUE, na.rm = TRUE, empty.cats = TRUE)
merge_temp1 = mergeProfilesLists(temp1, temp2, profNames = c("Gastric", "Ventricle"))
merge_temp2 = mergeProfilesLists(temp1, temp3, profNames = c("Ventricle", "Lung"))

plotProfiles(c(temp1, temp2, temp3), percentage=FALSE, multiplePlots = FALSE, HORIZVERT = TRUE,
             legendText = TRUE)
























