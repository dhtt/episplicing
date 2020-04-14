library('stringr', quietly=TRUE)
library("data.table", quietly=TRUE)
library("dplyr", quietly=TRUE)
library("VennDiagram", quietly=TRUE)

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
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()

E003_E005_res.sig <- as.data.frame(E003_E005_res[E003_E005_res$padj < 0.05 & !is.na(E003_E005_res$padj),])
E003_E005_res.sig_genes <- E003_E005_res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()

E004_E005_res.sig <- as.data.frame(E004_E005_res[E004_E005_res$padj < 0.05 & !is.na(E004_E005_res$padj),])
E004_E005_res.sig_genes <- E004_E005_res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()

all_genes = unique(E004_E005_res$groupID)

venn.plot <- venn.diagram(x = list(
  "E003_E004" = E003_E004_res.sig_genes$groupID,
  "E003_E005" = E003_E005_res.sig_genes$groupID,
  "E004_E005" = E004_E005_res.sig_genes$groupID
), 
filename=NULL, col="black", fill = c("blue", "red", "purple"), margin=0.05, alpha = 0.5)
jpeg("E003E004E005.jpg"); 
grid.draw(venn.plot); 
dev.off()
