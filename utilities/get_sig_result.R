library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(stringr)
library(boot)
library(stats)
library(parallel)
library(NMF)
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

#===========================RESULT LOADING + PREP ====================================
#Load files
all_res.exp = readRDS("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/report/all_DEU_res.exp.RDS")
all_res.his = readRDS("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/report/all_DEU_res.his.RDS")
all_res.met = readRDS("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/report/all_DEU_res.met.RDS")

combine_all_res <- function(all_res.exp, all_res.his, all_res.met){
  all_res = vector("list", 8)
  for (i in 1:length(all_res.his)){
    all_res[[i]] = all_res.his[[i]][seq(1, length(all_res.his[[i]]), 2)]
  }
  all_res[[7]] = all_res.met[seq(1, length(all_res.met), 2)]
  all_res[[8]] = all_res.exp[seq(1, length(all_res.exp), 2)]
  return(all_res)
}
all_res = combine_all_res(all_res.exp, all_res.his, all_res.met)
length(all_res) #all_res: 1-6 his, 7 met, 8 exp
lapply(all_res, length)

combine_all_res_exon <- function(all_res.exp, all_res.his, all_res.met){
  all_res = vector("list", 8)
  for (i in 1:length(all_res.his)){
    all_res[[i]] = all_res.his[[i]][seq(2, length(all_res.his[[i]]), 2)]
  }
  all_res[[7]] = all_res.met[seq(2, length(all_res.met), 2)]
  all_res[[8]] = all_res.exp[seq(2, length(all_res.exp), 2)]
  return(all_res)
}
all_res_exon = combine_all_res_exon(all_res.exp, all_res.his, all_res.met)
length(all_res_exon) #all_res: 1-6 his, 7 met, 8 exp
lapply(all_res_exon, length)


#Tissue spec 
get_tissue_spec_ref <- function(){
  RPKM = read.csv("~/Documents/BIOINFO/Episplicing/files/flank/57epigenomes.RPKM.pc", row.names=1, sep="")
  head(RPKM)
  epigenomes = c("E003", "E004", "E005", "E006", "E007", "E011","E012","E013", "E016", "E024", "E053","E054", "E065","E066","E071","E079","E094","E095", "E096", "E098", "E100","E105","E106", "E109","E113") #E022-E027
  RPKM = RPKM[, epigenomes]
  norm_RPKM = transpose(as.data.frame(apply(RPKM, 1, function(x) (x/max(x)))))
  head(norm_RPKM)
  TSI_RPKM = as.data.frame(apply(norm_RPKM, 1, function(x)((length(norm_RPKM)-sum(as.numeric(x)))/length(norm_RPKM))))
  TSI_RPKM$gene = rownames(RPKM)
  head(TSI_RPKM)
  ENS_list = TSI_RPKM[TSI_RPKM[[1]] >= 0.75, 2]
  
  ENS_gene_list = read.delim("~/Documents/BIOINFO/Episplicing/files/flank/Ensembl_v65.Gencode_v10.ENSG.gene_info.txt", header=FALSE)
  head(ENS_gene_list)
  colnames(ENS_gene_list) = c("ens", "chr", "start", "end", "strand", "feature", "symbol", "name")
  TSI_symbols = ENS_gene_list[ENS_gene_list$ens %in% ENS_list, "symbol"]
  return(TSI_symbols)
}
TSI_symbols = get_tissue_spec_ref()
epigenomes_names = c("H1 Cells", "Mesendoderm", "Trophoblast",  "Mesenchyma", "Neuronal Progenitor Cells",  "Endoderm", "Ectoderm", "Mesoderm", "HUES64 Cells", "ES-UCSF4", "Cortex-derived Neurospheres", "Ganglion Eminence-derived Neurospheres",
                     "Aorta", "Liver", "Brain Hippocampus Middle", "Esophagus", "Gastric", "Left Ventricle", "Lung", "Pancreas", "Psoas Muscle",   "Right Ventricle", "Sigmoid Colon", "Small Intestine", "Spleen")
get_tissue_spec_array <- function(){
  epigenomes_types = list(
    c("E003", "E004", "E005", "E006", "E007", "E011","E012","E013", "E016", "E024", "E053","E054", "E065","E066","E071","E079","E094","E095", "E096", "E098", "E100","E105","E106", "E109","E113"),
    c("E003", "E004", "E005", "E006", "E007", "E011","E012","E013", "E016", "E065","E066","E071","E079","E094","E095", "E096", "E098", "E100","E105","E106", "E109","E113"))
  tissue_type_list = vector("list")
  for (k in 1:length(epigenomes_types)){
    epigenomes = epigenomes_types[[k]]
    all_pairs = seq(1, length(epigenomes)*(length(epigenomes)-1)/2)
    print(all_pairs)
    idx_list = vector("list")
    idx = 1
    counter = 1
    n_epigenomes = length(epigenomes)
    while (counter <= n_epigenomes){
      print("new")
      print(counter)
      if (counter != n_epigenomes){
        idx_list[[counter]] = all_pairs[idx : (idx + (n_epigenomes - counter) - 1)]
      }
      else { idx_list[[counter]] = c(n_epigenomes*(n_epigenomes - 1)/2)}
      
      idx = idx + (n_epigenomes - counter)
      if (counter != 1){
        for (i in 1:(counter - 1)){
          temp = idx_list[[i]][counter-i]
          idx_list[[counter]] = union(idx_list[[counter]], temp)
        }
      }
      counter = counter + 1
      print(idx_list[[counter - 1]] )
    }
    tissue_type_list[[k]] = idx_list
  }
  return(tissue_type_list)
}
tissue_type_list = get_tissue_spec_array()
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac", "Methylation", "Expression")


#=================== Combine for each tissue then filter by TSI ========================
#Get genes for each tissues
get_genes_for_tissue <- function(res_list){
  all_genes_joined = vector("list") #List of 6 histone, each has sig genes for 25 tissues
  for (i in 1:length(res_list)){ #for each histone type
    res = res_list[[i]]
    all_tissues = vector("list")
    if (i == 6){
      tissue_list = tissue_type_list[[2]] }
    else {
      tissue_list = tissue_type_list[[1]] }
    for (j in 1:length(tissue_list)){
      print(paste(i, j, sep= ', '))
      tissue_type = tissue_list[[j]]
      all_tissues[[j]] = Reduce(union, res[tissue_type])
    }
    all_genes_joined[[i]] = all_tissues
  }
  temp = all_genes_joined[[6]][10:22]
  all_genes_joined[[6]][10:13] = c("")
  all_genes_joined[[6]][13:25] = temp
  return(all_genes_joined)
}
all_genes_joined = get_genes_for_tissue(all_res)
length(all_genes_joined)
lapply(all_genes_joined, length)
lapply(all_genes_joined[[6]], length)

get_all_len_before <- function(all_genes_joined){
  all_len_before = vector("list")
  for (i in 1:length(all_genes_joined)){
    len = lapply(all_genes_joined[[i]], length) 
    print(length(len))
    all_len_before[[i]] = len
  }
  all_len_before[[6]][c(10,11,12)] = NaN
  all_len_before = as.data.frame(do.call(cbind, all_len_before))
  rownames(all_len_before) = epigenomes_names
  colnames(all_len_before) = histone_type_list
  return(all_len_before)
}
all_len_before = get_all_len_before(all_genes_joined)

get_genes_for_tissue_filtered <- function(all_genes_joined){
  for (i in 1:length(all_genes_joined)){
    print(paste("New histone type", histone_type_list[i], sep=' '))
    for (j in 1: length(all_genes_joined[[i]])){
      print(paste(i, j, sep=' '))
      gene_set = all_genes_joined[[i]][[j]]
      gene_set = gene_set[gene_set %in% TSI_symbols]
      all_genes_joined[[i]][[j]] = gene_set
      print(length(gene_set))
    }
  } #List of 6 histone, each has FILTERED sig genes for 25 tissues
  return(all_genes_joined)
}
all_genes_joined_filtered = get_genes_for_tissue_filtered(all_genes_joined)
length(all_genes_joined) #gene names
lapply(all_genes_joined, length)
lapply(all_genes_joined[[1]], length) 

get_all_len_after <- function(all_genes_joined){
  all_len_after = vector("list")
  for (i in 1:length(all_genes_joined)){
    len = lapply(all_genes_joined[[i]], length) 
    print(length(len))
    all_len_after[[i]] = len
  }
  all_len_after[[6]][c(10,11,12)] = NaN
  all_len_after = as.data.frame(do.call(cbind, all_len_after))
  rownames(all_len_after) = epigenomes_names
  colnames(all_len_after) = histone_type_list
  return(all_len_after)
}
all_len_after = get_all_len_after(all_genes_joined_filtered)




























#=================== Heat map of result===================================

temp = all_distances[[8]]

all_len = lapply(all_res.exp, length)
all_len = unlist(unname(all_len))
all_len[tissue_type_list[[1]][[1]]]

get_distance_matrix_allres <- function(all_res){
  all_dist_matrix = vector("list", length(all_res))
  for (j in 1:length(all_res)){
    print(j)
    res_list = all_res[[j]]
    all_len = lapply(res_list, length)
    all_len = unlist(unname(all_len))
    dist_matrix = vector("list")
    if (j != 6){
      for (i in 1:(length(tissue_type_list[[1]]) - 1)){
        row_val = c(rep(NA, i), all_len[tissue_type_list[[1]][[i]]][1:(25-i)] )
        dist_matrix[[i]] = row_val
      }
      dist_matrix = as.data.frame(do.call(rbind, dist_matrix))
      dist_matrix[25, ] = NA
      dist_matrix[lower.tri(dist_matrix)] = t(dist_matrix)[lower.tri(dist_matrix)]
      rownames(dist_matrix) = epigenomes_names
      colnames(dist_matrix) = epigenomes_names
      all_dist_matrix[[j]] = dist_matrix
    }
    else {
      for (i in 1:(length(tissue_type_list[[2]]) - 1)){
        row_val = c(rep(NA, i), all_len[tissue_type_list[[2]][[i]]][1:(22-i)] )
        dist_matrix[[i]] = row_val
      }
      dist_matrix = as.data.frame(do.call(rbind, dist_matrix))
      dist_matrix[22, ] = NA
      dist_matrix[lower.tri(dist_matrix)] = t(dist_matrix)[lower.tri(dist_matrix)]
      rownames(dist_matrix) = epigenomes_names[-c(10,11,12)]
      colnames(dist_matrix) = epigenomes_names[-c(10,11,12)]
      all_dist_matrix[[j]] = dist_matrix
    }
  }
  return(all_dist_matrix)
}
all_distances = get_distance_matrix_allres(all_res)
temp = all_distances[[6]]
max(all_distances[[1]], na.rm = TRUE)
min(unlist(lapply(all_distances, function(x) min(x, na.rm = TRUE))))


epigenomes_names
epigenomes_type = c("Primary culture", "ES Cell-dev", "ES Cell-dev", "ES Cell-dev", "ES Cell-dev", 
                    "ES Cell-dev",  "ES Cell-dev", "ES Cell-dev", "Primary culture", "Primary culture", 
                    "Primary culture", "Primary culture", "Primary tissue", "Primary tissue", "Primary tissue", 
                    "Primary tissue", "Primary tissue", "Primary tissue", "Primary tissue", "Primary tissue",
                    "Primary tissue", "Primary tissue", "Primary tissue", "Primary tissue", "Primary tissue")
epigenomes_group = c("ES Cell", "ES Cell-dev", "ES Cell-dev", "ES Cell-dev", "ES Cell-dev", 
                     "ES Cell-dev",  "ES Cell-dev", "ES Cell-dev", "ES Cell", "ES Cell", 
                     "Neurospheres", "Neurospheres", "Heart", "Other", "Brain", 
                     "Digestive", "Digestive", "Heart", "Other", "Other",
                     "Muscle", "Heart", "Digestive", "Digestive", "Other")
unique(epigenomes_type)
epigenomes_annot = data.frame(type = epigenomes_type, group = epigenomes_group, row.names = epigenomes_names)
epigenomes_colors = list(c("hotpink", "PaleTurquoise", "MediumPurple"),
                         c("pink", "plum","wheat", "salmon", "slateblue", "darkseagreen", "gold", "lightblue"))

#separate images
breaks = seq(0, 25000, 250)
color = 'RdYlBu2:101'
folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/utilities/figures/"
tiff(paste(folder, histone_type_list[[1]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[1]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[1]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[2]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[2]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[2]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[3]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[3]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[3]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[4]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[4]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[4]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[5]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[5]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[5]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[6]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[6]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main=paste(histone_type_list[[6]]), 
         annCol = epigenomes_annot[-c(10,11,12),], annRow = epigenomes_annot[-c(10,11,12),],
         annColors = list(epigenomes_colors[[1]], epigenomes_colors[[2]][-c(7)]), color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[7]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[7]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[7]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()
tiff(paste(folder, histone_type_list[[8]], ".tiff", sep=''), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[8]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[8]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, color = color, breaks = breaks)
dev.off()




#=================== Get tissue for figure 1 =================== 
#"FGFR2" in exp and H3K36me3
for (i in tissue_type_list[[1]][[5]]){ # second index is tissue of interest ie 15:brain, 4: mesenchyma
  gene_array = all_res[[5]][[i]]
  exon_array = all_res_exon[[5]][[i]]
  gene_array1 = all_res[[8]][[i]]
  exon_array1 = all_res_exon[[8]][[i]]
  print(i)
  print(exon_array[match("NCAM1", gene_array)])
  print(exon_array1[match("NCAM1", gene_array1)])
  # print(intersect(, exon_array1[match("NCAM1", gene_array1)]))
  # print(exon_array[match("NCAM1", gene_array)])
}
tissue_type_list[[1]]
epigenomes_names[[16]]

length(all_res[[8]][[71]])
length(all_res_exon[[8]][[71]])

temp = all_res[[8]][[71]]
temp ==  "FGFR2"
match("FGFR2", temp)
temp = data.frame(all_res[[8]][[71]], all_res_exon[[8]][[71]])


