library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(viridis)
library(ggplot2)
library(ggraph)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")
all_epigenes_chi_list = readRDS("tissue_spec_genes/all_epigenes_chi_list.RDS")
all_epigenes_spearman_list = readRDS("tissue_spec_genes/all_epigenes_cor_spearman_list.RDS")
all_epigenes_pearson_list = readRDS("tissue_spec_genes/all_epigenes_cor_pearson_list.RDS")

get_entrez_id <-function(gene_list){
  entrez_gene_list = select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}
prepare_ck <- function(gene_list, analysis, option="all"){
  gene_list = lapply(gene_list, get_entrez_id)
  names(gene_list) = c("All", "ESC", "Mes", "Tro", "Gas", "Ven", "Mus")
  print("Preparing BP cluster")
  ck_bp <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                          OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = 0.01,
                          pAdjustMethod = "fdr", readable =TRUE)
  print("Preparing MF cluster")
  ck_mf <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                          OrgDb='org.Hs.eg.db', ont = "MF", qvalueCutoff = 0.01,
                          pAdjustMethod = "fdr", readable =TRUE)
  print("Preparing CC cluster")
  ck_cc <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                          OrgDb='org.Hs.eg.db', ont = "CC", qvalueCutoff = 0.01,
                          pAdjustMethod = "fdr", readable =TRUE)
  saveRDS(ck_bp, paste(analysis, "_BP_GO.RDS", sep=''))
  saveRDS(ck_mf, paste(analysis, "_MF_GO.RDS", sep=''))
  saveRDS(ck_cc, paste(analysis, "_CC_GO.RDS", sep=''))
  return(list(ck_bp, ck_mf, ck_cc))
}

#==== Pearson ====
pearson_res = prepare_ck(all_epigenes_pearson_list, analysis = "pearson")
pearson_bp = pearson_res[[1]]
pearson_mf = pearson_res[[2]]
pearson_cc = pearson_res[[3]]
#==== Spearman ====
spearman_res = prepare_ck(all_epigenes_spearman_list, analysis = "spearman")
spearman_bp = spearman_res[[1]]
spearman_mf = spearman_res[[2]]
spearman_cc = spearman_res[[3]]

tiff("ck_mf.tiff", units="in", width=12, height = 14, res=300)
dotplot(spearman_mf, showCategory = 50) +
  scale_color_viridis(option = "D") +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("ck_bp.tiff", units="in", width=12, height = 14, res=300)
dotplot(spearman_bp, showCategory = 50) +
  scale_color_viridis(option = "D") +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("ck_cc.tiff", units="in", width=12, height = 14, res=300)
dotplot(spearman_cc, showCategory = 50) +
  scale_color_viridis(option = "D") +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt")) +
  guides(color = guide_legend(override.aes = list(fill = "white")))
dev.off()

#==== Chi====
chi_res = prepare_ck(all_epigenes_chi_list, analysis = "chi")
chi_bp = chi_res[[1]]
chi_mf = chi_res[[2]]
chi_cc = chi_res[[3]]

tiff("chi_mf.tiff", units="in", width=12, height = 14, res=300)
dotplot(chi_mf, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("chi_bp.tiff", units="in", width=12, height = 14, res=300)
dotplot(chi_bp, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("chi_cc.tiff", units="in", width=12, height = 14, res=300)
dotplot(chi_cc, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()
#======
tiff("map_chi_mf.tiff", units="in", width=14, height = 14, res=300)
emapplot(chi_mf, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

#================EXON===============
exon_res = prepare_ck(all_exon_list, analysis = "exon")
exon_bp = exon_res[[1]]
exon_mf = exon_res[[2]]
exon_cc = exon_res[[3]]

tiff("exon_mf.tiff", units="in", width=12, height = 14, res=300)
dotplot(exon_mf, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("exon_bp.tiff", units="in", width=12, height = 14, res=300)
dotplot(exon_bp, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("exon_cc.tiff", units="in", width=12, height = 14, res=300)
dotplot(exon_cc, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("map_exon_mf.tiff", units="in", width=14, height = 14, res=300)
emapplot(exon_mf, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

tiff("map_exon_bp.tiff", units="in", width=14, height = 14, res=300)
emapplot(exon_bp, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

tiff("map_exon_cc.tiff", units="in", width=14, height = 14, res=300)
emapplot(exon_cc, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()
