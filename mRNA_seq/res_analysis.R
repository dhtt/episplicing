library('stringr', quietly=TRUE)
library("data.table", quietly=TRUE)
library("dplyr", quietly=TRUE)
library("VennDiagram", quietly=TRUE)
library("ggplot2")
library("gtools")
library("viridis")
library("reshape")

#===== LOAD DEXSEQ RESULT =====
setwd("~/Documents/BIOINFO/Episplicing/files/Result/combine/expression")
E003_E004_res = read.csv("E003_E004_res.csv", header=TRUE, sep = "\t")
E003_E005_res = read.csv("E003_E005_res.csv", header=TRUE, sep = "\t")
E004_E005_res = read.csv("E004_E005_res.csv", header=TRUE, sep = "\t")

#===== SUMMARIZE GENES WITH EXON Padj < 0.05 =====
E003_E004_res.sig <- as.data.frame(E003_E004_res[E003_E004_res$padj < 0.05 & !is.na(E003_E004_res$padj),])
E003_E004_res.sig_genes <- E003_E004_res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ","), 
         n_sigexon = str_count(featureID_by_groupID, ',') + 1) %>%
  dplyr::select(-featureID) %>%
  unique()

E003_E005_res.sig <- as.data.frame(E003_E005_res[E003_E005_res$padj < 0.05 & !is.na(E003_E005_res$padj),])
E003_E005_res.sig_genes <- E003_E005_res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ","), 
         n_sigexon = str_count(featureID_by_groupID, ',') + 1) %>%
  dplyr::select(-featureID) %>%
  unique()

E004_E005_res.sig <- as.data.frame(E004_E005_res[E004_E005_res$padj < 0.05 & !is.na(E004_E005_res$padj),])
E004_E005_res.sig_genes <- E004_E005_res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ","), 
         n_sigexon = str_count(featureID_by_groupID, ',') + 1) %>%
  dplyr::select(-featureID) %>%
  unique()

all_genes = unique(E004_E005_res$groupID)
all_sig_genes = list(E003_E004_res.sig_genes, E003_E005_res.sig_genes, E004_E005_res.sig_genes)

table(E003_E004_res.sig_genes$n_sigexon >= 3)
table(E003_E005_res.sig_genes$n_sigexon >= 3)
table(E004_E005_res.sig_genes$n_sigexon >= 3)
#===== VISUALIZATION =====
##===== Venn diagram =====
venn.plot <- venn.diagram(x = list(
  "E003_E004" = E003_E004_res.sig_genes$groupID,
  "E003_E005" = E003_E005_res.sig_genes$groupID,
  "E004_E005" = E004_E005_res.sig_genes$groupID
), 
filename=NULL, col="black", fill = c("yellow", "red", "purple"), margin=0.05, alpha = 0.5)
jpeg("E003E004E005.jpg"); 
grid.draw(venn.plot); 
dev.off()

##===== Heatmap =====
epi_id1 = "E003"
epi_id2 = "E004"
epi_id3 = "E005"
pairs = transpose(as.data.frame(combn(c(epi_id1, epi_id2, epi_id3), 2)))
as_genes_id = as_genes$gene
cs_genes_id = constitutive_genes$gene


deu_genes_as = lapply(all_sig_genes, function(x) intersect(x[[1]], as_genes_id))
deu_genes_cs = lapply(all_sig_genes, function(x) intersect(x[[1]], cs_genes_id))
lapply(deu_genes_as, length)
lapply(deu_genes_cs, length)
as_gene_summary = cbind(pairs, AS_genes, length(as_genes_id)-AS_genes)
colnames(genes_summary) = c("epi_id1", "epi_id2", "AS_genes", "nonAS_genes")

rep_pairs = as.data.frame(cbind(c(epi_id1, epi_id2, epi_id3), c(epi_id1, epi_id2, epi_id3), 0, length(all_genes)))
colnames(rep_pairs) = colnames(genes_summary)
genes_summary = rbind(genes_summary, rep_pairs, .id = NULL)
temp = melt(genes_summary, id.vars = c("epi_id1","epi_id2"), variable_name = c("type"))

ggplot(data = temp, aes(x = epi_id1, y = epi_id2)) +
  geom_tile(aes(fill = value))  +
  facet_grid(. ~ type) +
  geom_tile(aes(fill = value))  +
  scale_fill_viridis(discrete=TRUE) +
  labs(title = "Number of alternatively spliced genes in pairwise comparison of epigenomes", x=NULL, y=NULL) +
  theme_bw()







