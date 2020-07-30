library(data.table, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(plyr)
library(rtracklayer)
library(ggplot2)
library(reshape2)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/flank")

#---------------------------------------------
folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/flank"
all_res_list.pearcor_p = readRDS(paste("all_res_list.pearcor_p.RDS", sep='/'))
all_res_list.pearcor_p.met = readRDS(paste("all_res_list.pearcor_p.met.RDS", sep='/'))
all_res_list.pearcor_p[[7]] = all_res_list.pearcor_p.met
lapply(all_res_list.pearcor_p, dim)

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
      else if (method == "fisher" | method == "chisq" | method == "pairedcor" | method == "pearcor_p" | method == "spearcor_p" ){
        all_res_sig[[j-1]] = all_res[all_res.col <= p_sig  & !is.na(all_res.col), 1]
      }
    }
    all_res_list_sig[[i]] = all_res_sig
  }
  return(all_res_list_sig)
}
all_res_list.pearcor_p_sig = get_all_res_list_sig(all_res_list.pearcor_p, "pearcor_p", p_sig=0.05)
lapply(all_res_list.pearcor_p_sig, length)
# ----2-----
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
histone_type_list = list("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac", "Methylation")

#START here
all_res_list.pearcor_p_sig_joinedtissue = vector("list") #List of 6 histone, each has sig genes for 25 tissues
for (i in 1:length(all_res_list.pearcor_p_sig)){ #for each histone type
  res = all_res_list.pearcor_p_sig[[i]]
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
  all_res_list.pearcor_p_sig_joinedtissue[[i]] = all_tissues
}
temp = all_res_list.pearcor_p_sig_joinedtissue[[6]][10:22]
all_res_list.pearcor_p_sig_joinedtissue[[6]][10:13] = c("")
all_res_list.pearcor_p_sig_joinedtissue[[6]][13:25] = temp
length(all_res_list.pearcor_p_sig_joinedtissue)
lapply(all_res_list.pearcor_p_sig_joinedtissue, length)
lapply(all_res_list.pearcor_p_sig_joinedtissue[[7]], length)

all_len_before = vector("list")
for (i in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
  len = lapply(all_res_list.pearcor_p_sig_joinedtissue[[i]], length) 
  print(length(len))
  all_len_before[[i]] = len
}
all_len_before[[6]][c(10,11,12)] = NaN
all_len_before = as.data.frame(do.call(cbind, all_len_before))
epigenomes_names = c("H1 Cells", "Mesendoderm", "Trophoblast",  "Mesenchyma", "Neuronal Progenitor Cells",  "Endoderm", "Ectoderm", "Mesoderm", "HUES64 Cells", "ES-UCSF4", "Cortex-derived Neurospheres", "Ganglion Eminence-derived Neurospheres",
                     "Aorta", "Liver", "Brain Hippocampus Middle", "Esophagus", "Gastric", "Left Ventricle", "Lung", "Pancreas", "Psoas Muscle",   "Right Ventricle", "Sigmoid Colon", "Small Intestine", "Spleen")
rownames(all_len_before) = epigenomes_names
colnames(all_len_before) = histone_type_list

for (i in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
  print(paste("New histone type", histone_type_list[i], sep=' '))
  for (j in 1: length(all_res_list.pearcor_p_sig_joinedtissue[[i]])){
    print(paste(i, j, sep=' '))
    gene_set = all_res_list.pearcor_p_sig_joinedtissue[[i]][[j]]
    gene_set = gene_set[gene_set %in% TSI_symbols]
    all_res_list.pearcor_p_sig_joinedtissue[[i]][[j]] = gene_set
    print(length(gene_set))
  }
} #List of 6 histone, each has FILTERED sig genes for 25 tissues
length(all_res_list.pearcor_p_sig_joinedtissue) #gene names
lapply(all_res_list.pearcor_p_sig_joinedtissue, length)
lapply(all_res_list.pearcor_p_sig_joinedtissue[[1]], length) 

all_len_after = vector("list")
for (i in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
  len = lapply(all_res_list.pearcor_p_sig_joinedtissue[[i]], length) 
  print(length(len))
  all_len_after[[i]] = len
}
all_len_after[[6]][c(10,11,12)] = NaN
all_len_after = as.data.frame(do.call(cbind, all_len_after))
rownames(all_len_after) = epigenomes_names
colnames(all_len_after) = histone_type_list

#END here
length(all_res_list.pearcor_p_sig_joinedtissue[[1]])
paste(all_res_list.pearcor_p_sig_joinedtissue[[4]][[12]], collapse = ', ')
temp = Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[7]])
paste(temp, collapse = ', ')
temp = union(all_res_list.pearcor_p_sig_joinedtissue[[7]][[12]], all_res_list.pearcor_p_sig_joinedtissue[[2]][[12]])

# -----3-----
length(intersect(all_res_list.pearcor_p_sig_joinedtissue[[1]][[1]], all_res_list.pearcor_p_sig_joinedtissue[[1]][[2]]))
Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]])

all_distances = vector("list")
for (i in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
  temp = all_res_list.pearcor_p_sig_joinedtissue[[i]]
  distances = sapply(temp, function(x) sapply(temp, function(y) length(intersect(x, y))))
  diag(distances) = NA
  distances = distances/max(distances, na.rm = TRUE)
  rownames(distances) = epigenomes_names
  colnames(distances) = epigenomes_names
  if (i == 6){
    distances = distances[-c(10, 11, 12), -c(10, 11, 12)]
  }
  all_distances[[i]] = distances
}

library(NMF)
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
tiff("pairwise_distance2.tiff", width = 30, height = 40, units="in", res=100)
par(mfrow = c(4,2))
for (i in 1:7){
  if (i != 6){
    aheatmap(all_distances[[i]], Rowv = TRUE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
             annColors = epigenomes_colors)
  }
  else {
    aheatmap(all_distances[[i]], Rowv = TRUE, Colv = TRUE, scale="none", 
             main=paste(histone_type_list[[i]]), 
             annCol = epigenomes_annot[-c(10,11,12),], annRow = epigenomes_annot[-c(10,11,12),],
             annColors = list(epigenomes_colors[[1]], epigenomes_colors[[2]][-c(7)]))
  }
}
dev.off()

#separate images
folder = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/utilities/figures"
breaks = c(seq(0, 0.2, 0.2), seq(0.21, 1, 0.01))
color = '-RdYlBu2:82'
tiff(paste(folder, histone_type_list[[1]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[1]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[1]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[2]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[2]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[2]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[3]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[3]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[3]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[4]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[4]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[4]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[5]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[5]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[5]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[6]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[6]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main=paste(histone_type_list[[6]]), 
         annCol = epigenomes_annot[-c(10,11,12),], annRow = epigenomes_annot[-c(10,11,12),],
         annColors = list(epigenomes_colors[[1]], epigenomes_colors[[2]][-c(7)]),
         breaks = breaks, color = color)
dev.off()
tiff(paste(folder, histone_type_list[[7]], sep='/'), width = 13, height = 11, units="in", res=100)
aheatmap(all_distances[[7]], Rowv = TRUE, Colv = TRUE, scale="none", 
         main = paste(histone_type_list[[7]]), annCol = epigenomes_annot, annRow = epigenomes_annot,
         annColors = epigenomes_colors, breaks = breaks, color = color)
dev.off()
# -----4-----
epigenomes_names[c(2,8,7,3)] #H1 1
epigenomes_names[c(9, 11, 12)] #H1 2
epigenomes_names[c(11, 12, 9, 10)] #H2
epigenomes_names[c(14,13,18)] #H2 2
epigenomes_names[c(16, 20, 25, 19)] #H3 1
epigenomes_names[c(11, 12, 9, 10)] #H3 2
epigenomes_names[c(11, 12, 9, 10,1)] #H4 1
epigenomes_names[c(3, 4, 6, 2)] #H4 2
epigenomes_names[c(16,23)]
epigenomes_names[c(8,3,7,6,4)] #H5
epigenomes_names[c(14,25,18)] #H6 1
epigenomes_names[c(24,23,16,19)] #H6 1
epigenomes_names[c(20, 22, 21)] #H6 1
paste0(Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]][c(8,3,7,6,4)]), collapse = "|")

Reduce(intersect, all_res_list.pearcor_p_sig_joinedtissue[[5]])

# -----3-----
all_len_after.his = all_len_after[,1:6]
all_len_after.his = as.data.frame(lapply(all_len_after.his, function(x) as.numeric(as.character(x))))
rownames(all_len_after.his) = epigenomes_names
colnames(all_len_after.his) = histone_type_list[1:6]
head(all_len_after.his)

n_tissues = 25
all_res_list_sig_joined = vector("list")
for (i in 1:n_tissues){
  all_res_sig = vector("list")
  for (j in 1:length(all_res_list.pearcor_p_sig_joinedtissue)){
    print( paste(j, i, sep = ", "))
    all_res_sig[[j]] = all_res_list.pearcor_p_sig_joinedtissue[[j]][[i]]
  }
  temp = Reduce(union, all_res_sig)
  all_res_list_sig_joined[[i]] = temp
}
length(all_res_list_sig_joined)
lapply(all_res_list_sig_joined, length)
lapply(all_res_list.pearcor_p_sig_joinedtissue[[6]][[10]], length)
all_res_list_sig_joined[[1]]

all_len_after.his$`Total genes with DHP (overlap)` = rowSums(all_len_after.his, na.rm=TRUE)
all_len_after.his$`Total genes with DHP (non-overlap)` = as.numeric(lapply(all_res_list_sig_joined, length))
sd_nan <- function(x){
  return(sd(x, na.rm = TRUE))
}
temp = all_len_after.his %>% 
  summarise_all(sd_nan) 

all_len_after.his$`Total genes with DMP (overlap)` = all_len_before$Methylation
all_len_after.his$`Total genes with DMP (non-overlap)` = all_len_after$Methylation
# -----4-----
library(formattable)
tissue_formatter <- formatter("span", 
                              style = x ~ style(
                                width = suffix(x, "px"),
                                font.weight = "bold", 
                                color = ifelse(x == "Total", "black", "gray")))
formattable(all_len_after.his,
            # align =c("l", "l", "c","c","c","c"), 
            list(
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K4me3` = color_tile("white", "wheat"),
              `H3K9me3` = color_tile("white", "wheat"),
              `H3K27me3` = color_tile("white", "wheat"),
              `H3K36me3` = color_tile("white", "wheat"), 
              `H3K4me1` = color_tile("white", "wheat"),
              `H3K27ac` = color_tile("white", "wheat"),
              `Total genes with DHP (overlap)` = color_tile("white", "LightSalmon"),
              `Total genes with DHP (non-overlap)` = color_tile("white", "LightSalmon"),
              `Total genes with DMP (overlap)` = color_tile("white", "LightSalmon"),
              `Total genes with DMP (non-overlap)` = color_tile("white", "LightSalmon")
            )
)


