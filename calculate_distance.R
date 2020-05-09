library(data.table, quietly=TRUE)
library(stringr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(GenomicRanges)

#Get DEU middle
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine2")
epi_id = c("E003","E004","E005","E094","E095","E096")

#Load reference genome with promoter and exon
refgen = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/mrna_seq/reference_genome.whole.sorted.gtf"
refgen = import.gff(refgen)
refgen = as(refgen, "GRanges")
refgen = refgen[refgen$type == "exonic_part" | refgen$type == "promoter",]
refgen$exonic_part_number[is.na(refgen$exonic_part_number)] = "000"
refgen$exonic_part_number = paste('E', refgen$exonic_part_number, sep='')
refgen$middle = round(abs(start(refgen)-end(refgen))/2 + start(refgen), digits = 0)
refgen.df = as.data.frame(refgen)
head(refgen.df)
getwd()
#Load DEXSeq result with DEU exon
get_pos_DEU_exon <-function(epi_id1, epi_id2, refgen.df){
  file_path = normalizePath(paste(getwd(), "expression", "annotatedcounts", sep='/'))
  pair = paste(epi_id1, epi_id2, sep='_')
  dexseq_res = list.files(file_path, pattern=pair, full.names=TRUE)[[1]]
  dexseq_res = read.csv(dexseq_res, header=TRUE, sep = "\t")
  colnames(dexseq_res) = c("gene", "exon", "exonBaseMean", "stat", "pvalue", "padj", "log2fold")
  dexseq_res <- as.data.frame(dexseq_res[dexseq_res$padj < 0.05 & !is.na(dexseq_res$padj)
                                         & !is.na(dexseq_res$log2fold),])
  
  dexseq_res.sig = as.vector(c(paste(dexseq_res$gene, dexseq_res$exon , sep =':'),
                               paste(unique(dexseq_res$gene), "E000", sep=':')))
  dexseq_res.siggenes = unique(dexseq_res$gene)
  
  exp.sig_pos = refgen.df %>%
    filter(gene_id %in% sig_gene) %>%
    mutate(gene = gene_id, 
           ID = paste(gene_id, exonic_part_number , sep =':'),
           range = paste(start, end, sep=',')) %>%
    filter(ID %in% dexseq_res.sig) %>%
    dplyr::select(gene, exonic_part_number, ID, range, middle) %>%
    group_by(gene) %>%
    mutate(middle = paste(middle, collapse = ","),
           range = paste(range, collapse = ",")) %>%
    dplyr::select(gene, range, middle) %>%
    unique()
  return(exp.sig_pos)
}
pair1_2.deu_pos = get_pos_DEU_exon(epi_id[1], epi_id[2], refgen.df)
pair1_3.deu_pos = get_pos_DEU_exon(epi_id[1], epi_id[3], refgen.df)
pair2_3.deu_pos = get_pos_DEU_exon(epi_id[2], epi_id[3], refgen.df)
pair4_5.deu_pos = get_pos_DEU_exon(epi_id[4], epi_id[5], refgen.df)
pair4_6.deu_pos = get_pos_DEU_exon(epi_id[4], epi_id[6], refgen.df)
pair5_6.deu_pos = get_pos_DEU_exon(epi_id[5], epi_id[6], refgen.df)
pair_list.exp_pos = list(pair1_2.deu_pos, pair1_3.deu_pos, pair2_3.deu_pos,
                         pair4_5.deu_pos, pair4_6.deu_pos, pair5_6.deu_pos)


#Load MAnorm peaks
get_hpeak_pos <- function(type, epi_id1, epi_id2){
  print("---> Prepare file")
  file_path = normalizePath(paste(getwd(), type, "annotatedcounts", sep='/'))
  pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
               paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
  hisdiff.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  hisdiff.exon = fread(hisdiff.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff.pi = fread(hisdiff.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff_tables_list = list(hisdiff.exon, hisdiff.pi)
  
  print("---> Processing")
  for (i in range(1:2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table <- hisdiff_table %>%
      # filter(grepl(gene, V9, fixed = TRUE)) %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
             V15 = if_else(is.na(V15), 0, V15), V16 = if_else(is.na(V16), 0, V16),
             V13 = as.numeric(as.character(V13)), V14 = as.numeric(as.character(V14)),
             V13 = if_else(is.na(V13), 0, V13), V14 = if_else(is.na(V14), 1, V14),
             V17 = V9) %>%
      dplyr::select(c(1,3,7,11,12,13,14,15,16,17))
    colnames(hisdiff_table) = c("chr", "feature", "strand", "start", "end","m-val","p-val", epi_id1, epi_id2, "gene")
    if (i == 1){
      hisdiff_table$feature = "exon"
    }
    hisdiff_tables_list[[i]] = hisdiff_table
  }
  
  print("---> Finding middle pos")
  hisdiff_tables_list = rbind(hisdiff_tables_list[[1]], hisdiff_tables_list[[2]])
  colnames(hisdiff_tables_list) = c("chr", "feature", "strand", "start", "end","m-val","p-val", epi_id1, epi_id2, "gene")
  hisdiff_tables_list$gene = as.data.table(str_split_fixed(hisdiff_tables_list$gene, "\\;|\\ ",8))[[2]]
  hisdiff_tables_list$middle = round(abs(hisdiff_tables_list$start - hisdiff_tables_list$end)/2 + hisdiff_tables_list$start, digits = 0)
  
  print("---> Summarizing peaks middle pos")
  hpeak.middle = hisdiff_tables_list %>% 
    mutate(gene = gsub('\"', '', gene)) %>%
    filter(`p-val` <= 0.05, abs(`m-val`) >= 1, gene %in% sig_gene) %>%
    group_by(gene) %>%
    mutate(middle = paste(middle, collapse = ',')) %>%
    dplyr::select(gene, middle) %>%
    unique()
  return(hpeak.middle)
}

# pair1_2.his_pos = get_hpeak_pos("H3K36me3", epi_id[1], epi_id[2])
# pair1_3.his_pos = get_hpeak_pos("H3K36me3", epi_id[1], epi_id[3])
# pair2_3.his_pos = get_hpeak_pos("H3K36me3", epi_id[2], epi_id[3])
# pair4_5.his_pos = get_hpeak_pos("H3K36me3", epi_id[4], epi_id[5])
# pair4_6.his_pos = get_hpeak_pos("H3K36me3", epi_id[4], epi_id[6])
# pair5_6.his_pos = get_hpeak_pos("H3K36me3", epi_id[5], epi_id[6])

pair1_2.his_pos = get_hpeak_pos("H3K27ac", epi_id[1], epi_id[2])
pair1_3.his_pos = get_hpeak_pos("H3K27ac", epi_id[1], epi_id[3])
pair2_3.his_pos = get_hpeak_pos("H3K27ac", epi_id[2], epi_id[3])
pair4_5.his_pos = get_hpeak_pos("H3K27ac", epi_id[4], epi_id[5])
pair4_6.his_pos = get_hpeak_pos("H3K27ac", epi_id[4], epi_id[6])
pair5_6.his_pos = get_hpeak_pos("H3K27ac", epi_id[5], epi_id[6])
pair_list.his_pos = list(pair1_2.his_pos, pair1_3.his_pos, pair2_3.his_pos,
                         pair4_5.his_pos, pair4_6.his_pos, pair5_6.his_pos)

get_all_pair_distance <- function(pair_list.exp_pos, pair_list.his_pos){
  all_pair_distance = vector("list", length(pair_list.exp_pos))
  for (i in 1:length(pair_list.exp_pos)){
    print(paste("New pair ", i, sep=''))
    pair.deu_pos = pair_list.exp_pos[[i]]
    pair.his_pos = pair_list.his_pos[[i]]
    common_genes= intersect(pair.deu_pos$gene, pair.his_pos$gene)
    pair_distance = vector("list", length(common_genes))
    for (j in 1:length(common_genes)){
      gene = common_genes[[j]]
      exp.range = as.numeric(strsplit(pair.deu_pos[pair.deu_pos$gene == gene,"range"][[1]], ',')[[1]])
      his.middle = as.numeric(strsplit(pair.his_pos[pair.his_pos$gene == gene,"middle"][[1]], ',')[[1]])
      
      pair_distance[[j]] = get_min_dis(exp.range, his.middle)
    }
    pair_distance = Reduce(c, pair_distance)
    all_pair_distance[[i]] = pair_distance
  }
  all_pair_distance = data.frame(Reduce(c, all_pair_distance))
  colnames(all_pair_distance) = c("distance")
  return(all_pair_distance)
}
all_pair_distance_H3K36me3 = get_all_pair_distance(pair_list.exp_pos, pair_list.his_pos)
all_pair_distance_H3K27ac
all_pair_distance_H3K36me3$type = "H3K36me3"
all_pair_distance_H3K27ac$type = "H3K27ac"

all_pair_distance = rbind(all_pair_distance_H3K36me3, all_pair_distance_H3K27ac)
table(all_pair_distance$distance > 100000)
histogram.color = "#39038f"

tiff("d_histone.tiff", units="in", width=8, height = 5, res=300)
ggplot(all_pair_distance, aes(x = distance, fill = type, color = type)) +
  geom_histogram(bins = 200, alpha = 0.7) +
  scale_x_continuous(breaks = round(seq(0, max(all_distance.collapse$distance), by = 2000), 0),
                     limits = c(0, 25000),
                     name = "Distance (bp)") +
  scale_y_continuous(breaks = round(seq(0, 2000, by = 200), 0),
                     limits = c(0, 1600),
                     name = "Counts") +
  scale_color_viridis_d(option = "D") +
  scale_fill_viridis_d(option = "D") +
  ggtitle("Distance between significant histone peaks to nearest promoter/exons") +
  theme_minimal() +
  theme(
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text.x=element_text(colour="black", size = 9),
        axis.text.y=element_text(colour="black", size = 9),
        plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()


