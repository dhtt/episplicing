start_time <- Sys.time()
#===== LOAD PACKAGES ======
library('stringr', quietly=TRUE)
library("data.table", quietly=TRUE)
library("dplyr", quietly=TRUE)
library("DEXSeq", quietly=TRUE)
library("optparse", quietly=TRUE)
suppressPackageStartupMessages( library( "DEXSeq" ) )

option_list = list(
  make_option(c("-f", "--countfolder"), type="character", default="~/Documents/BIOINFO/Episplicing/files/Result/combine/expression", 
              help="path to folder of counts", metavar="character"),
  make_option(c("-a", "--epigenome1"), type="character", default=NULL, 
              help="ID of first epigenome", metavar="character"),
  make_option(c("-b", "--epigenome2"), type="character", default=NULL, 
              help="ID of second epigenome", metavar="character"),
  make_option(c("-g", "--referencegenome"), type="character", default="/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/mrna_seq/reference_genome.gtf", 
              help="path to flattened reference genome", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#===== PREPARE DATA =====
setwd(opt$countfolder)
inDir = normalizePath(getwd())
epi_id1 = opt$epigenome1
epi_id2 = opt$epigenome2
pair = paste(paste('^', epi_id1, ".*count.txt$", sep=''), paste('^',  epi_id2, ".*count.txt$", sep=''), sep='|')
count_files = list.files(inDir, pattern=pair, full.names=TRUE)
file_names = as.data.table(str_split_fixed(basename(count_files), "\\_", 3))
gtf_files = opt$referencegenome

print(paste("---> Working folder: ", opt$countfolder, sep=''))
print("---> Count files: ")
print(basename(count_files))
print(paste("---> Reference genome: ", gtf_files, sep=''))

sampleTable = data.frame(
  row.names = c(file_names$V2),
  condition = c(file_names$V1))

#===== RUN DEXSEQ =====
print("---> Inputting to DEXSeq")
dxd = DEXSeqDataSetFromHTSeq(
  count_files,
  sampleData = sampleTable,
  design = ~ sample + exon + condition:exon,
  flattenedfile= normalizePath(gtf_files)
)

print("---> Getting DEXSeq result")
dxd.res = DEXSeq(dxd, quiet = FALSE)

#===== SAVING RESULTS =====
print("---> Saving DEXSeq normalized counts")
dxd.count = data.frame(counts(dxd.res, normalized = TRUE))
colnames(dxd.count) = paste(file_names$V1, file_names$V2, sep='_')
normedcount_name = paste(paste(epi_id1, epi_id2, sep='_'), "normedcount.csv", sep='_')
write.table(dxd.count, normedcount_name, quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
# dxd.count = read.csv("temp_count.csv", header=TRUE, sep = ",")

print("---> Saving DEXSeq result")
result_name = paste(paste(epi_id1, epi_id2, sep='_'), "res.csv", sep='_')
write.table(as.data.frame(dxd.res[c(1,2,3,5,6,7,10)]), result_name, 
            quote=FALSE, sep="\t", dec=".", row.names=FALSE, col.names=TRUE)
# dxd.res = read.csv(result_name, header=TRUE, sep = ",")

print("===> FINISHED!")
end_time <- Sys.time()
end_time - start_time

