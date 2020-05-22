library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(data.table)
library(dplyr)
background.title = "#39038f"
col.title = "white"
col.highlight = "#fffa70"
fontsize.title = 14
type.color = c("#EE442F", "#63ACBE")
margin = c(20, 20)

gen = "hg19"
#======================================================================================
example_gene_folder = normalizePath("/Users/dhthutrang/example_E005_E100_ORAI3")
gene_name = "ORAI3"
epi_id1 = "E005"
epi_id2 = "E100"
epi_name1 = "Trophoblast"
epi_name2 = "Right Ventricle"
all_files = list.files(example_gene_folder, full.names = TRUE)

#===== EXP =====
transcript_file = all_files[grep("NCBI", all_files)]
transcript = import.gff(transcript_file)
transcript = as(transcript, "GRanges")
transcript = transcript[transcript$type == "exon",]
chr <- as.character(unique(seqnames(transcript)))[[1]]
itrack <- IdeogramTrack(genome = gen, chromosome = chr, name = gene_name)
gtrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(transcript,
                           genome = gen, chromosome = chr,
                           name = "Transcripts", 
                           transcriptAnnotation = "transcript",
                           fontsize.title = fontsize.title,
                           col = "black", fill="#85c0f9"
)

transcript_sig_file = all_files[grep("res.csv", all_files)]
transcript_sig = fread(transcript_sig_file)
transcript_sig = transcript_sig[transcript_sig$padj <= 0.05, c(6,11,12)]
colnames(transcript_sig) = c("padj", "start", "end")
# transcript_sig = transcript_sig[c(2,3,4,5,6,9,10),]
if (dim(transcript_sig)[[1]] != 0){
  grtrack_sig <- HighlightTrack(trackList = grtrack, chromosome = chr,
                                start = transcript_sig$start - 0,
                                end = transcript_sig$end + 0,
                                col = "magenta", fill = "magenta", alpha = 0.3,
                                inBackground = FALSE)
} else {
  grtrack_sig = grtrack
}

plotTracks(list(itrack, gtrack, grtrack, grtrack_sig), extend.left = 2500, extend.right = 2500)


#===== MET =====
prepare_metdiff_example <- function(epi_id1, epi_id2, gene){
  metdiff.exon = fread(normalizePath(all_files[grep('.*exon.*normedratio.csv.txt', all_files, fixed = FALSE)]), 
                       header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff.pi = fread(normalizePath(all_files[grep('.*pi.*normedratio.csv.txt', all_files, fixed = FALSE)]), 
                     header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  head(metdiff.exon)
  head(metdiff.pi)
  metdiff_tables_list = list(metdiff.exon, metdiff.pi)

  for (i in range(1:2)){
    metdiff_table = metdiff_tables_list[[i]]
    metdiff_table <- metdiff_table %>%
      filter(V3 == "exonic_part" | V3 == "promoter" | V3 == "intron") %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V14 = as.numeric(as.character(V14)), V15 = as.numeric(as.character(V15)),
             V14 = if_else(is.na(V14), 0, V14), V15 = if_else(is.na(V15), 0, V15),
             V16 = gene) %>%
      dplyr::select(c(1,3,7,11,12,14,15,16))

    colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", epi_id1, epi_id2, "gene")
    metdiff_tables_list[[i]] = metdiff_table
  }
  metdiff_tables_list = rbind(metdiff_tables_list[[1]], metdiff_tables_list[[2]])
  metdiff_tables_list = as(metdiff_tables_list, "GRanges")
  return(metdiff_tables_list)
}
met_count = prepare_metdiff_example(epi_id1, epi_id2, gene_name)
mTrack <- DataTrack(met_count, name = "Methylation",
                    type = c("a"),
                    col = type.color,
                    legend = TRUE)

get_sig_met <- function(epi_id1, epi_id2, gene, thres_metdiff){
  metdiff = fread(normalizePath(all_files[grep('.*diff.txt.txt', all_files, fixed = FALSE)]), 
                       header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  metdiff_table <- metdiff %>%
        mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
               V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
               V15 = if_else(is.na(V15), 1, V15), V16 = if_else(is.na(V16), 0, V16),
               V17 = gene) %>%
        dplyr::select(c(1,3,7,11,12,15,16,17))
  colnames(metdiff_table) = c("chr", "feature", "strand", "start", "end", "q-val", "metdiff","gene")
  metdiff_table = metdiff_table[metdiff_table$`q-val` <= 0.05 & abs(metdiff_table$`metdiff`) >= thres_metdiff, ]
  metdiff_table = as(metdiff_table, "GRanges")
  metdiff_table = as.data.frame(metdiff_table)
  return(metdiff_table)
}
met_sig = get_sig_met(epi_id1, epi_id2, gene_name, 70)

if (dim(met_sig)[[1]] != 0){
  mTrack_sig <- HighlightTrack(trackList = mTrack, chromosome = chr,
                               start = met_sig$start - 100,
                               end = met_sig$end + 200,
                               col = "#fffa70", fill = "#fffa70")
} else {
  mTrack_sig = mTrack
}

plotTracks(list(itrack, gtrack, grtrack_sig, mTrack_sig),
           extend.left = 2500, extend.right = 2500,
           groups = c(epi_id1, epi_id2),
           background.title = "#597ca8")

#===== HIS =====
his_type = c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3")
his_type = c("H3K36me3", "H3K9me3", "H3K4me3", "H3K27me3", "H3K27ac", "H3K4me1")
prepare_hisdiff_example <- function(epi_id1, epi_id2, gene, type){
  file_id_exon = paste(".*exon", paste(type, "txt", sep='.'), sep='.*')
  file_id_pi = paste(".*pi", paste(type, "txt", sep='.'), sep='.*')
  hisdiff.exon = fread(normalizePath(all_files[grep(file_id_exon, all_files, fixed = FALSE)]), 
                       header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff.pi = fread(normalizePath(all_files[grep(file_id_pi, all_files, fixed = FALSE)]), 
                     header=FALSE, stringsAsFactors = FALSE, quote=FALSE)
  hisdiff_tables_list = list(hisdiff.exon, hisdiff.pi)

  for (i in range(1:2)){
    hisdiff_table = hisdiff_tables_list[[i]]
    hisdiff_table <- hisdiff_table %>%
      filter(V3 == "exonic_part" | V3 == "promoter" | V3 == "intron") %>%
      mutate(V11 = if_else(V11 == -1, V4, V11), V12 = if_else(V12 == -1, V5, V12),
             V15 = as.numeric(as.character(V15)), V16 = as.numeric(as.character(V16)),
             V15 = if_else(is.na(V15), 0, V15), V16 = if_else(is.na(V16), 0, V16),
             V13 = as.numeric(as.character(V13)), V14 = as.numeric(as.character(V14)),
             V13 = if_else(is.na(V13), 0, V13), V14 = if_else(is.na(V14), 1, V14),
             V17 = V9) %>%
      dplyr::select(c(1,3,7,11,12,13,14,15,16,17))

    colnames(hisdiff_table) = c("chr", "feature", "strand", "start", "end","m-val","p_val", epi_id1, epi_id2, "gene")
    hisdiff_tables_list[[i]] = hisdiff_table
  }
  hisdiff_tables_list = rbind(hisdiff_tables_list[[1]], hisdiff_tables_list[[2]])
  hisdiff_tables_list = as(hisdiff_tables_list, "GRanges")
  return(hisdiff_tables_list)
}
get_all_his_tables <- function(){
  all_his_counts = vector("list", length(his_type))
  all_his_sigs = vector("list", length(his_type))
  for (i in 1:length(his_type)){
    type = his_type[i]
    his_table = prepare_hisdiff_example(epi_id1, epi_id2, gene_name, type)
    his_table.count = his_table[, c(epi_id1, epi_id2)]
    his_table.sig = as.data.frame((his_table[abs(his_table$`m-val`) >= 1 & his_table$`p_val` <= 0.05]))
    all_his_counts[[i]] = his_table.count
    all_his_sigs[[i]] = his_table.sig
  }
  return(list(all_his_counts, all_his_sigs))
}
all_his_counts = get_all_his_tables()
all_his_sigs = all_his_counts[[2]]
all_his_counts = all_his_counts[[1]]
lapply(all_his_sigs, function(x) dim(x)[[1]])

get_all_his_track.sigs <- function(){
  all_his_track.sigs = vector("list", length(all_his_counts))
  for (i in 1:length(all_his_counts)){
    his_count = all_his_counts[[i]]
    his_sig = all_his_sigs[[i]]
    if (dim(his_sig)[[1]] != 0){
      hTrack <- DataTrack(his_count, name = his_type[i],
                          type = c("a"), legend = FALSE,
                          col = type.color)
      hTrack.sig <- HighlightTrack(trackList = hTrack, chromosome = chr,
                                   start = his_sig$start  - 100,
                                   end = his_sig$end  + 100,
                                   col = col.highlight, fill = col.highlight)
      all_his_track.sigs[[i]] = hTrack.sig
    }
    else {
      hTrack <- DataTrack(his_count, name = his_type[i],
                          type = c("a"), legend = FALSE,
                          col = type.color)
      all_his_track.sigs[[i]] = hTrack
    }
  }
  return(all_his_track.sigs)
}
all_his_track.sigs = get_all_his_track.sigs()

temp = all_his_sigs[[6]]
sum(temp$end - temp$start)/(max(end(transcript)) - min(start(transcript)))*100

#===== ALL TRACKS =====
all_tracks = list(itrack, gtrack, grtrack_sig)
all_tracks = append(all_tracks, all_his_track.sigs)
all_tracks = append(all_tracks, mTrack_sig)
plotTracks(all_tracks,
           extend.left = 2500, extend.right = 2500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           fontsize = fontsize.title,
           col.title = col.title,
           groups = c(epi_name1, epi_name2),
           margin = margin, 
           stackHeight = 0.5)

tiff(paste(paste(epi_id1, epi_id2, gene_name, sep='_'), "tiff", sep='.'), 
     units="in", width=9, height=10, res=300)
plotTracks(all_tracks,
           extend.left = 500, extend.right = 500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           fontsize = fontsize.title,
           col.title = col.title,
           groups = c(epi_name1, epi_name2),
           margin = margin, 
           stackHeight = 0.5)
dev.off()





