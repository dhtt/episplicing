library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(data.table)
library(dplyr)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine2")

SYNPO.gtf = import.gff("SYNPO.gtf") #from NCBI
SYNPO.gtf = as(SYNPO.gtf, "GRanges")
SYNPO.gtf = SYNPO.gtf[SYNPO.gtf$type == "exon",]

gtrack <- GenomeAxisTrack()
gen = "hg19"
chr <- as.character(unique(seqnames(SYNPO.gtf)))
itrack <- IdeogramTrack(genome = gen, chromosome = chr, name = "SYNPO")

head(SYNPO.gtf)

# background.panel = ""
background.title = "#39038f"
col.title = "white"
col.highlight = "#fffa70"
fontsize.title = 11
type.color = c("#EE442F", "#63ACBE")
gene_name = "SYNPO"

grtrack <- GeneRegionTrack(SYNPO.gtf, 
                           genome = gen, chromosome = chr,
                           name = "Transcripts",
                           transcriptAnnotation = "transcript",
                           col = "black", fill="#85c0f9"
                           )
plotTracks(list(itrack, gtrack, grtrack),
           extend.left = 2500, extend.right = 2500)

#=====Read his=====
prepare_hisdiff_example <- function(type, epi_id1, epi_id2, gene){
  gene = paste(gene, '\"', sep = '')
  file_path = normalizePath(paste(getwd(), type, "annotatedcounts", sep='/'))
  pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
               paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
  hisdiff.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  hisdiff.exon = fread(hisdiff.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff.pi = fread(hisdiff.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  print(head(hisdiff.pi$V9))
  hisdiff_tables_list = list(hisdiff.exon, hisdiff.pi)
  
  for (i in range(1:2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table <- hisdiff_table %>%
      filter(grepl(gene, V9, fixed = TRUE)) %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
             V15 = if_else(is.na(V15), 0, V15), V16 = if_else(is.na(V16), 0, V16),
             V13 = as.numeric(as.character(V13)), V14 = as.numeric(as.character(V14)),
             V13 = if_else(is.na(V13), 0, V13), V14 = if_else(is.na(V14), 1, V14),
             V17 = V9) %>%
      dplyr::select(c(1,3,7,11,12,13,14,15,16,17))
    
    colnames(hisdiff_table) = c("chr", "feature", "strand", "start", "end","m-val","p_val", epi_id1, epi_id2, "gene")
    if (i == 1){
      hisdiff_table$feature = "exon"
    }
    hisdiff_tables_list[[i]] = hisdiff_table
  }
  return(hisdiff_tables_list)
}
SYNPO.H3K27ac = prepare_hisdiff_example("H3K27ac", "E094", "E095", 'SYNPO')
SYNPO.H3K27ac = rbind(SYNPO.H3K27ac[[1]], SYNPO.H3K27ac[[2]])
SYNPO.H3K27ac = as(SYNPO.H3K27ac, "GRanges")
SYNPO.H3K27ac.counts = SYNPO.H3K27ac[, c("E094","E095")]
#ADJUST THRESHOLD OFR HISTONE DIFERENTLY -> top N REGIONS(here n = 14)
SYNPO.H3K27ac.sig = as.data.frame((SYNPO.H3K27ac[abs(SYNPO.H3K27ac$`m-val`) >= 2 & SYNPO.H3K27ac$`p_val` <= 0.05]))


SYNPO.H3K36me3 = prepare_hisdiff_example("H3K36me3", "E094", "E095", "SYNPO")
SYNPO.H3K36me3 = rbind(SYNPO.H3K36me3[[1]], SYNPO.H3K36me3[[2]])
SYNPO.H3K36me3 = as(SYNPO.H3K36me3, "GRanges")
SYNPO.H3K36me3.counts = SYNPO.H3K36me3[, c("E094","E095")]
#ADJUST THRESHOLD OFR HISTONE DIFERENTLY -> top N REGIONS(here n = 14)
SYNPO.H3K36me3.sig = as.data.frame((SYNPO.H3K36me3[SYNPO.H3K36me3$`p_val` <= 0.05]))

hTrack1 <- DataTrack(SYNPO.H3K27ac.counts, name = "H3K27ac", 
                     type = c("a"), legend = FALSE,
                     col = type.color)
hTrack2 <- DataTrack(SYNPO.H3K36me3.counts, name = "H3K36me3", 
                     type = c("a"), legend = FALSE,
                     col = type.color)

ht1 <- HighlightTrack(trackList = hTrack1, chromosome = 5, 
                      start = SYNPO.H3K27ac.sig$start  - 1000, 
                      width = SYNPO.H3K27ac.sig$width  + 1000,
                      col = col.highlight, fill = col.highlight)
ht2 <- HighlightTrack(trackList = hTrack2, chromosome = 5, 
                      start = SYNPO.H3K36me3.sig$start - 1000, 
                      width = SYNPO.H3K36me3.sig$width + 1000,
                      col = col.highlight, fill = col.highlight)
# ht1 <- HighlightTrack(trackList = dTrack1, chromosome = 5,
#                       start = c(112231319, 112220770),
#                       end=c(112232684, 112221422))
# ht2 <- HighlightTrack(trackList = dTrack2, chromosome = 5,
#                       start = c(112213796),
#                       end = c(112214829))

plotTracks(list(itrack, gtrack, grtrack, ht1, ht2),
           extend.left = 2500, extend.right = 2500,
           groups = c("E094","E095"))
#======Read met  =========
prepare_metdiff_example <- function(epi_id1, epi_id2, gene){
  gene = paste(gene, '\"', sep = '')
  file_path = normalizePath(paste(getwd(), "methylation", "annotatedcounts", sep='/'))
  pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
               paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
  metdiff.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  metdiff.exon = fread(metdiff.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff.pi = fread(metdiff.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_tables_list = list(metdiff.exon, metdiff.pi)
  
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table <- metdiff_table %>%
      filter(grepl(gene, V9)) %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V14 = as.numeric(as.character(V14)), V15 = as.numeric(as.character(V15)),
             V14 = if_else(is.na(V14), 0, V14), V15 = if_else(is.na(V15), 0, V15),
             V16 = gene) %>%
      dplyr::select(c(1,3,7,11,12,14,15,16))
    
    colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", epi_id1, epi_id2, "gene")
    if (i == 1){
      metdiff_table$feature = "exon"
    }
    metdiff_tables_list[[i]] = metdiff_table
  }
  return(metdiff_tables_list)
}
get_sig_met <- function(epi_id1, epi_id2, gene, thres_metdiff){
  gene = paste(gene, '\"', sep = '')
  file_path = normalizePath(paste(getwd(), "methylation", "annotateddiff", sep='/'))
  pair = paste(paste('^pi', paste(epi_id1, epi_id2, sep='_'), sep='_'),
               paste('^exon', paste(epi_id1, epi_id2, sep='_'), sep='_'), sep='|')
  metdiff.diff = list.files(file_path, pattern=pair, full.names=TRUE)
  metdiff.exon = fread(metdiff.diff[[1]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff.pi = fread(metdiff.diff[[2]], header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_tables_list = list(metdiff.exon, metdiff.pi)
  
  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table <- metdiff_table %>%
      filter(grepl(gene, V9)) %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
             V15 = if_else(is.na(V15), 1, V15), V16 = if_else(is.na(V16), 0, V16),
             V17 = gene) %>%
      dplyr::select(c(1,3,7,11,12,15,16,17))
    
    colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", "q-val", "metdiff","gene")
    if (i == 1){
      metdiff_table$feature = "exon"
    }
    metdiff_table = metdiff_table[metdiff_table$`q-val` <= 0.05 & abs(metdiff_table$`metdiff`) >= thres_metdiff, ]
    metdiff_tables_list[[i]] = metdiff_table
  }
  return(metdiff_tables_list)
}

SYNPO.met = prepare_metdiff_example("E094", "E095", "SYNPO")
SYNPO.met = rbind(SYNPO.met[[1]], SYNPO.met[[2]])
SYNPO.met = as(SYNPO.met, "GRanges")
mTrack <- DataTrack(SYNPO.met, name = "Methylation", 
                    type = c("smooth"), span = 0.02,
                    col = type.color,
                    legend = TRUE)

SYNPO.met.sig = get_sig_met("E094","E095", "SYNPO", 50)
SYNPO.met.sig = rbind(SYNPO.met.sig[[1]], SYNPO.met.sig[[2]])
SYNPO.met.sig = as(SYNPO.met.sig, "GRanges")
SYNPO.met.sig = as.data.frame(SYNPO.met.sig)

ht3 <- HighlightTrack(trackList = mTrack, chromosome = 5, 
                      start = SYNPO.met.sig$start - 1000, 
                      width = SYNPO.met.sig$width + 1000,
                      col = "#fffa70", fill = "#fffa70")
plotTracks(list(grtrack, ht1, ht2, ht3, gtrack, itrack),
           # extend.left = 2500, extend.right = 2500,
           groups = c("E094","E095"),
           background.title = "#597ca8")

#======Read exp========
SYNPO.exp = fread("SYNPO_exp.gtf")
SYNPO.exp = SYNPO.exp[,c("seqnames","start", "end", "width", "strand", "E094", "E095", "padj")]
colnames(SYNPO.exp) = c("chr","start", "end", "width", "strand", "E094", "E095", "padj")
SYNPO.exp$padj = as.numeric(as.character(SYNPO.exp$padj))
SYNPO.exp$gene = "SYNPO"
# SYNPO.exp$start = SYNPO.exp$start + SYNPO.exp$width/2
SYNPO.exp = as(SYNPO.exp, "GRanges")
SYNPO.exp.counts = SYNPO.exp[, c("E094","E095")]
eTrack <- DataTrack(SYNPO.exp.counts, name = "Expression",
                    type = c("h"), 
                    fill = "85c0f9")

SYNPO.exp.sig = SYNPO.exp[, "padj"]
SYNPO.exp.sig = SYNPO.exp[SYNPO.exp$padj <= 0.05 & !is.na(SYNPO.exp$padj)]
sig_range4 = as.data.frame(SYNPO.exp.sig)
ht4 <- HighlightTrack(trackList = grtrack, chromosome = 5,
                      start = sig_range4$start - 50,
                      width = sig_range4$width + 50,
                      col = "orange", fill = "#fffa70")
ht4 <- HighlightTrack(trackList = grtrack, chromosome = 5,
                      start = c(150029866, 150031904, 150035966) - 0,
                      end = c(150031903, 150033794, 150038769) + 0,
                      col = "magenta", fill = "magenta", alpha = 0.3,
                      inBackground = FALSE)

tiff("SYNPO.tiff", units="in", width=8, height= 5, res=300)
plotTracks(list(itrack, gtrack,ht4, ht1, ht2, ht3),
           extend.left = 2500, extend.right = 2500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           col.title = col.title,
           groups = c("E094","E095"),
           margin = c(20, 20))
dev.off()
# temp = as.data.frame(SYNPO.exp.counts)
# 
# temp$start + temp$width/2
