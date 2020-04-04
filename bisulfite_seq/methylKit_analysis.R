#===== LOAD PACKAGES =====
library("genomation", quietly=TRUE)
library("data.table", quietly=TRUE)
library("methylKit", quietly=TRUE)
library("optparse", quietly=TRUE)

option_list = list(
  make_option(c("-f", "--countfolder"), type="character", default="/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/Result/combine/methylation", 
              help="path to folder of counts", metavar="character"),
  make_option(c("-a", "--epigenome1"), type="character", default=NULL, 
              help="ID of first epigenome", metavar="character"),
  make_option(c("-b", "--epigenome2"), type="character", default=NULL, 
              help="ID of second epigenome", metavar="character"),
  make_option(c("-c", "--mincov"), type="integer", default=10, 
              help="minimum coverage to be included", metavar="character"),
  make_option(c("-n", "--numcores"), type="integer", default=1, 
              help="number of processing cores", metavar="character")
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

print(paste("---> Working folder: ", opt$countfolder, sep=''))
print("---> Count files: ")
print(basename(count_files))

mincov = opt$mincov
cores = opt$numcores
#===== RUN METHYLKIT =====
#TODO: CORES
methylrawlist = methRead(as.list(count_files),
                         sample.id = list(epi_id1, epi_id2), treatment = c(0,1), mincov=mincov,
                         assembly = 'hg19', header=FALSE, context="CpG", resolution="base",
                         pipeline=list(fraction=TRUE, chr.col=1, start.col=2, end.col=2, 
                                       strand.col=3, freqC.col=4, coverage.col=5))

normalized_methylrawlist = normalizeCoverage(methylrawlist, method = "median", adjust = "fdr") 
normalized_methylrawlist.united = unite(normalized_methylrawlist, destrand=FALSE, mc.cores = cores)

#===== SAVING RESULTS =====
print("---> Saving MethylKit result")
diffmeth = calculateDiffMeth(normalized_methylrawlist.united, adjust='fdr', mc.cores = cores, 
                             save.db = TRUE, suffix = paste(epi_id1, epi_id2, sep='_'))
# diffmeth.df = data.frame(diffmeth)

print("---> Saving MethylKit normalized methylation ratio")
ratio1 = normalized_methylrawlist.united$numCs1/normalized_methylrawlist.united$coverage1
ratio2 = normalized_methylrawlist.united$numCs2/normalized_methylrawlist.united$coverage2
normedratio = data.frame(cbind(getData(normalized_methylrawlist.united)[,1:4], ratio1, ratio2))
colnames(normedratio) = c(colnames(normedratio)[1:4], c(epi_id1, epi_id2))
normedratio_name = paste(paste(epi_id1, epi_id2, sep='_'), "normedratio.csv", sep='_')
write.table(normedratio, normedratio_name, quote=FALSE, sep=",", dec=".", row.names=TRUE, col.names=TRUE)
# normedratio = read.csv("E003_E004_normedratio.csv", header=TRUE, sep = ",")

print("===> FINISHED!")