library(RRreg)
library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(psych)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")

#---------------------------------------------
all_res_list.pearcor = readRDS("all_res_con/all_res_list.pearcor.RDS")
all_res_list.pearcor_p = readRDS("all_res_con/all_res_list.pearcor_p.RDS")
all_res_list.spearcor = readRDS("all_res_con/all_res_list.spearcor.RDS")
all_res_list.bootcor_fisher = readRDS("all_res_con/all_res_list.bootcor_fisher.RDS")
# all_res_list.chisq = analyze_array_list(all_pairs.exp, all_pairs.his_list, "chisq")
# all_res_list.randcor = analyze_array_list(all_pairs.exp, all_pairs.his_list, "randcor")

# ------------ Get significant results ------------
# ----1-----
get_all_res_list_sig <- function(all_res_list, method, r_sig=0.5, p_sig= 0.05){
  all_res_list_sig = vector("list", length(all_res_list))
  for (i in 1:length(all_res_list)) {
    all_res = all_res_list[[i]]
    all_res_sig = vector("list", ncol(all_res)-1)
    for (j in 2:ncol(all_res)){
      all_res.col = as.numeric(all_res[[j]])
      if (method == "randcor" | method == "pearcor" | method == "spearcor"){
        all_res_sig[[j-1]] = all_res[abs(all_res.col) >= r_sig  & !is.na(all_res.col), 1] 
      }
      else if (method == "bootcor_fisher"){
        all_res_sig[[j-1]] = all_res[abs(all_res.col) >= r_sig & abs(all_res.col) <= 1 & !is.na(all_res.col), 1] 
      }
      else if (method == "fisher" | method == "chisq" | method == "pairedcor" | method == "pearcor_p" ){
        all_res_sig[[j-1]] = all_res[all_res.col <= p_sig  & !is.na(all_res.col), 1]
      }
    }
    all_res_list_sig[[i]] = all_res_sig
  }
  return(all_res_list_sig)
}

all_res_list.pearcor_sig = get_all_res_list_sig(all_res_list.pearcor, "pearcor", r_sig=0.5)
all_res_list.pearcor_p_sig = get_all_res_list_sig(all_res_list.pearcor_p, "pearcor_p", p_sig=0.1)
all_res_list.spearcor_sig = get_all_res_list_sig(all_res_list.spearcor, "spearcor", r_sig=0.5)
all_res_list.bootcor_fisher_sig = get_all_res_list_sig(all_res_list.bootcor_fisher, "bootcor_fisher", r_sig=0.5)
# all_res_list.chisq_sig = get_all_res_list_sig(all_res_list.chisq, "chisq", p_sig = 0.1)
# all_res_list.randcor_sig = get_all_res_list_sig(all_res_list.randcor, "randcor", r_sig=0.7)

# ----2-----
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
n_pairs = 15
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27ac", "H3K27me3", "H3K36me3")

all_res_list.pearcor_sig_joined = join_all_res_list_sig(all_res_list.pearcor_sig)
all_res_list.pearcor_p_sig_joined = join_all_res_list_sig(all_res_list.pearcor_p_sig)
all_res_list.spearcor_sig_joined = join_all_res_list_sig(all_res_list.spearcor_sig)
all_res_list.bootcor_fisher_sig_joined = join_all_res_list_sig(all_res_list.bootcor_fisher_sig)
# all_res_list.chisq_sig_joined = join_all_res_list_sig(all_res_list.chisq_sig)
# all_res_list.randcor_sig_joined = join_all_res_list_sig(all_res_list.randcor_sig)

# ----3-----
get_tissue_specific_gene <- function(all_epigenes_list, pair_list, method="intersect"){
  if (method == "union"){
    all_epigenes_intersect = Reduce(union, all_epigenes_list[pair_list])
  }
  else {
    all_epigenes_intersect = Reduce(intersect, all_epigenes_list[pair_list])
  }
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
# saveRDS(all_res_list.pearcor_sig_joined_genes, "all_res_con/tissue_spec/all_res_list.pearcor_sig_joined_genes.RDS")
all_res_list.pearcor_p_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.pearcor_p_sig_joined, x))
# saveRDS(all_res_list.pearcor_p_sig_joined_genes, "all_res_con/tissue_spec/all_res_list.pearcor_p_sig_joined_genes.RDS")
all_res_list.spearcor_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.spearcor_sig_joined, x))
# saveRDS(all_res_list.spearcor_sig_joined_genes, "all_res_con/tissue_spec/all_res_list.spearcor_sig_joined_genes.RDS")
all_res_list.bootcor_fisher_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.bootcor_fisher_sig_joined, x))
# saveRDS(all_res_list.bootcor_fisher_sig_joined_genes, "all_res_con/tissue_spec/all_res_list.bootcor_fisher_sig_joined_genes.RDS")

lapply(all_res_list.pearcor_p_sig_joined_genes, function(x) length(x))
lapply(all_res_list.pearcor_p_sig_joined_genes, function(x) paste(x, collapse = ', '))
temp = lapply(all_res_list.spearcor_sig_joined_genes, function(x) setdiff(x, all_res_list.pearcor_sig_joined_genes[[1]]))
lapply(temp, function(x) paste(x, collapse = ', '))


#======== PAIRED COR =======
gene_freq = as.data.table(fread("exp_id.txt", header=FALSE))
gene_freq = gene_freq %>%
  group_by(V1) %>%
  summarise(n = n())
gene_freq = gene_freq$n
all_genes = fread("gene_id.txt", header = FALSE)
compare_paired_r <- function(all_res_list.pearcor, all_res_list.bootcor_fisher){
  all_res_list.pairedcor_sig = vector("list", length(all_res_list.pearcor))
  for (i in 1:length(all_res_list.pearcor)){
    dcor1 = as.matrix(all_res_list.pearcor[[i]])
    dcor2 = as.matrix(all_res_list.bootcor_fisher[[i]])
    dcor = matrix(NA, ncol=ncol(dcor1)-1, nrow=nrow(dcor2))
    
    for(j in 1:nrow(dcor1)){
      for(k in 2:ncol(dcor1)){
        n_cor = as.numeric(gene_freq[j])
        if(!is.na(dcor1[j,k]) & !is.na(dcor2[j,k]) & n_cor > 3){
          dcor[j,k-1] <- paired.r(as.numeric(dcor1[j,k]), as.numeric(dcor2[j,k]), n = n_cor)$p
        }
        else {
          dcor[j,k-1] <- NA
        }
      }
    }
    dcor = as.data.table(dcor)
    dcor = cbind(all_genes$V1, dcor)
    all_res_list.pairedcor_sig[[i]] = as.data.frame(dcor)
  }
  return(all_res_list.pairedcor_sig)
}

all_res_list.pairedcor = compare_paired_r(all_res_list.pearcor, all_res_list.bootcor_fisher)
saveRDS(all_res_list.pairedcor, "all_res_con/all_res_list.pairedcor.RDS")
all_res_list.pairedcor_sig = get_all_res_list_sig(all_res_list.pairedcor, "pairedcor", p_sig = 0.1)
all_res_list.pairedcor_sig_joined = join_all_res_list_sig(all_res_list.pairedcor_sig)
all_res_list.pairedcor_sig_joined_genes = lapply(tissue_type_list, function(x) get_tissue_specific_gene(all_res_list.pairedcor_sig_joined, x, method = "union"))
# saveRDS(all_res_list.pairedcor_sig_joined_genes, "all_res_con/tissue_spec/all_res_list.pairedcor_sig_joined_genes.RDS")

lapply(all_res_list.pairedcor_sig_joined_genes, function(x) length(x))
lapply(all_res_list.pairedcor_sig_joined_genes, function(x) paste(x, collapse = ', '))

##===== COMBINE =====
all_genes_pearcor = Reduce(union, all_res_list.pearcor_sig_joined)
all_genes_pearcor_p = Reduce(union, all_res_list.pearcor_p_sig_joined)
all_genes_spearcor = Reduce(union, all_res_list.spearcor_sig_joined_genes)
all_genes_bootcor_fisher = Reduce(union, all_res_list.bootcor_fisher_sig_joined)
all_genes_pairedcor = Reduce(union, all_res_list.pairedcor_sig_joined)

library(VennDiagram)
venn.plot <- venn.diagram(x = list(
  "Pearson cor" = all_genes_pearcor,
  "Pearson-p cor" = all_genes_pearcor_p,
  "Spearman cor" = all_genes_spearcor,
  # "Bootstrapped cor " = all_genes_bootcor,
  "Bootstrapped-Fisher cor " = all_genes_bootcor_fisher,
  "Paired cor " = all_genes_pairedcor
  # "Fisher" = all_genes_fisher
), 
filename=NULL, col="black", fill = c("orange", "red", "purple", "pink", "blue"), margin=0.1, alpha = 0.4,
cat.dist = 0.1, print.mode = "raw")
jpeg("compare_cat5.jpg", res = 300, width = 7, height = 7, units="in"); 
grid.draw(venn.plot); 
dev.off()

