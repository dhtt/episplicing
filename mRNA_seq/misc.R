#Venn diagram for DEX
# dxd_stranded_E003E004.res
# dxd_stranded_E003E005.res
# dxd_stranded_E004E005.res

# dxd_stranded_E003E004.res_ordered = dxd_stranded_E003E004.res[order(dxd_stranded_E003E004.res$padj),]
dxd_stranded_E003E004.res.sig <- as.data.frame(dxd_stranded_E003E004.res[dxd_stranded_E003E004.res$padj < 0.05 & !is.na(dxd_stranded_E003E004.res$padj),])
dxd_stranded_E003E004.res.genes <- dxd_stranded_E003E004.res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()
#-----
# dxd_stranded_E003E005.res_order = dxd_stranded_E003E005.res[order(dxd_stranded_E003E005.res$padj),]
dxd_stranded_E003E005.res.sig <- as.data.frame(dxd_stranded_E003E005.res[dxd_stranded_E003E005.res$padj < 0.05 & !is.na(dxd_stranded_E003E005.res$padj),])
dxd_stranded_E003E005.res.genes <- dxd_stranded_E003E005.res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()
#-----
# dxd_stranded_E004E005.res_ordered = dxd_stranded_E004E005.res[order(dxd_stranded_E004E005.res$padj),]

dxd_stranded_E004E005.res.sig <- as.data.frame(dxd_stranded_E004E005.res[dxd_stranded_E004E005.res$padj < 0.05 & !is.na(dxd_stranded_E004E005.res$padj),])
dxd_stranded_E004E005.res.genes <- dxd_stranded_E004E005.res.sig %>%
  dplyr::select(groupID, featureID) %>%
  group_by(groupID) %>%
  mutate(featureID_by_groupID = paste(featureID, collapse = ",")) %>%
  dplyr::select(-featureID) %>%
  unique()

venn.plot <- venn.diagram(x = list(
  "E003_E004" = dxd_stranded_E003E004.res.genes$groupID,
  "E003_E005" = dxd_stranded_E003E005.res.genes$groupID,
  "E004_E005" = dxd_stranded_E004E005.res.genes$groupID
), 
filename=NULL, col="black", fill = c("blue", "red", "purple"), margin=0.05, alpha = 0.5, print.mode="percent")
jpeg("venn_jpeg.jpg"); 
grid.draw(venn.plot); 
dev.off();