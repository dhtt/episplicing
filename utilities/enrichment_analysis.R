library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(viridis)
library(ggplot2)
library(ggraph)
setwd("/Users/dhthutrang/Documents/BIOINFO/Episplicing/files/gene_level")
get_entrez_id <-function(gene_list){
  entrez_gene_list = select(org.Hs.eg.db,
                            keys = gene_list,
                            columns = c("ENTREZID", "SYMBOL"),
                            keytype = "SYMBOL")
  return(entrez_gene_list[[2]])
}
prepare_ck <- function(gene_list, option, q_val = 0.05){
  gene_list = lapply(gene_list, get_entrez_id)
  # names(gene_list) = c("ESC", "Mes", "Tro", "Gas", "Ven", "Mus")
  names(gene_list) = c("All", "ESC", "Mes", "Tro", "Gas", "Ven", "Mus")
  if (option == "BP"){
    print("Preparing BP cluster")
    ck_bp <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                            OrgDb='org.Hs.eg.db', ont = "BP", qvalueCutoff = q_val,
                            pAdjustMethod = "fdr", readable =TRUE)
    return(ck_bp)
  }
  else if (option == "MF"){
    print("Preparing MF cluster")
    ck_mf <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                            OrgDb='org.Hs.eg.db', ont = "MF", qvalueCutoff = q_val,
                            pAdjustMethod = "fdr", readable =TRUE)
    return(ck_mf)
  }
  else if (option == "CC"){
    print("Preparing CC cluster")
    ck_cc <- compareCluster(geneCluster = gene_list, fun = "enrichGO",
                            OrgDb='org.Hs.eg.db', ont = "CC", qvalueCutoff = q_val,
                            pAdjustMethod = "fdr", readable =TRUE)
    return(ck_cc)
  }
}
make_dotplot <- function(res, annot_type, type){
  dplot = dotplot(res, showCategory = 50) +
    scale_color_viridis(option = "D", alpha = 0.8) +
    ggtitle(paste("Enriched GO Terms for Tissue-specific Gene Sets", annot_type, type, sep='_')) +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      axis.text.x=element_text(colour="black", size = 9),
      axis.text.y=element_text(colour="black", size = 9),
      plot.margin = unit(c(20,20,20,20), "pt"))
  return(dplot)
}
make_emapplot <- function(res, type, annot_type){
  eplot = emapplot(res, pie="count", pie_scale=1.5, layout="kk", showCategory = 20, legend_n=5)+
    ggtitle(paste("Enriched GO Terms for Tissue-specific Gene Sets", annot_type, type, sep='_')) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
      plot.margin = unit(c(50, 50, 50, 50), "pt")) +
    scale_fill_manual(drop=FALSE,values = plasma(7))
  return(eplot)
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

tiff("ck_mf.tiff", units="in", width=12, height = 20, res=200)
dotplot(spearman_mf, showCategory = 50) +
  scale_color_viridis(option = "D") +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("ck_bp.tiff", units="in", width=12, height = 20, res=200)
dotplot(spearman_bp, showCategory = 50) +
  scale_color_viridis(option = "D") +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("ck_cc.tiff", units="in", width=12, height = 20, res=200)
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

tiff("chi_mf.tiff", units="in", width=12, height = 20, res=200)
dotplot(chi_mf, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("chi_bp.tiff", units="in", width=12, height = 20, res=200)
dotplot(chi_bp, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("chi_cc.tiff", units="in", width=12, height = 20, res=200)
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
tiff("map_chi_mf.tiff", units="in", width=20, height = 20, res=200)
emapplot(chi_mf, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

#================EXON===============
exon_res = prepare_ck(all_exon_list, analysis = "exon")
exon_bp = exon_res[[1]]
exon_mf = exon_res[[2]]
exon_cc = exon_res[[3]]

tiff("exon_mf.tiff", units="in", width=12, height = 20, res=200)
dotplot(exon_mf, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("exon_bp.tiff", units="in", width=12, height = 20, res=200)
dotplot(exon_bp, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("exon_cc.tiff", units="in", width=12, height = 20, res=200)
dotplot(exon_cc, showCategory = 50) +
  scale_color_viridis(option = "D", alpha = 0.8) +
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.text.x=element_text(colour="black", size = 9),
    axis.text.y=element_text(colour="black", size = 9),
    plot.margin = unit(c(20,20,20,20), "pt"))
dev.off()

tiff("map_exon_mf.tiff", units="in", width=20, height = 20, res=200)
emapplot(exon_mf, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Molecular Function)") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

tiff("map_exon_bp.tiff", units="in", width=20, height = 20, res=200)
emapplot(exon_bp, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Biological Process)") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

tiff("map_exon_cc.tiff", units="in", width=20, height = 20, res=200)
emapplot(exon_cc, pie="count", pie_scale=1.5, layout="nicely", showCategory = 20, legend_n=5)+
  ggtitle("Enriched GO Terms for Tissue-specific Gene Sets (Cellular Component)") +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.margin = unit(c(50, 50, 50, 50), "pt")) +
  scale_fill_manual(drop=FALSE,values = plasma(7))
dev.off()

#=================================== CONTINUOUS SCALE ===================================
folder = "all_res_con"
pearcor_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.pearcor_sig_joined_genes.RDS", sep='/'))
spearcor_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.spearcor_sig_joined_genes.RDS", sep='/'))
pearcor_p_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.pearcor_p_sig_joined_genes.RDS", sep='/'))
spearcor_p_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.spearcor_p_sig_joined_genes.RDS", sep='/'))
# chisq_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.chisq_sig_joined_genes.RDS", sep='/'))
bootcor_fisher_res_list = readRDS(paste(folder, "tissue_spec/all_res_list.bootcor_fisher_sig_joined_genes.RDS", sep='/'))
# pairedcor_res_list = readRDS("all_res_con/tissue_spec/all_res_list.pairedcor_sig_joined_genes.RDS")

pearcor_res.mf = prepare_ck(pearcor_res_list, option = "MF")
spearcor_res.mf = prepare_ck(spearcor_res_list, option = "MF")
pearcor_p_res.mf = prepare_ck(pearcor_p_res_list, option = "MF")
spearcor_p_res.mf = prepare_ck(spearcor_p_res_list, option = "MF")
# chisq_res.mf = prepare_ck(chisq_res_list, option = "MF")
bootcor_fisher_res.mf = prepare_ck(bootcor_fisher_res_list, option = "MF")
# pairedcor_res.mf = prepare_ck(pairedcor_res_list, option = "MF")
all_res.mf = list(pearcor_res.mf, spearcor_res.mf, pearcor_p_res.mf, spearcor_p_res.mf, bootcor_fisher_res.mf)

pearcor_res.bp = prepare_ck(pearcor_res_list, option = "BP")
spearcor_res.bp = prepare_ck(spearcor_res_list, option = "BP")
pearcor_p_res.bp = prepare_ck(pearcor_p_res_list, option = "BP")
spearcor_p_res.bp = prepare_ck(spearcor_p_res_list, option = "BP")
# chisq_res.bp = prepare_ck(chisq_res_list, option = "BP")
bootcor_fisher_res.bp = prepare_ck(bootcor_fisher_res_list, option = "BP")
# pairedcor_res.bp = prepare_ck(pairedcor_res_list, option = "BP")
all_res.bp = list(pearcor_res.bp, spearcor_res.bp, pearcor_p_res.bp, spearcor_p_res.bp, bootcor_fisher_res.mf)

cor_type = c("pearcor", "spearcor", "pearcor_p", "spearcor_p", "chisq")


#DOT PLOT
tiff(paste(folder, "annot_mf", "dotplot", paste(cor_type[1], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.mf[[1]], type = cor_type[1], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "dotplot", paste(cor_type[2], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.mf[[2]], type = cor_type[2], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "dotplot", paste(cor_type[3], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.mf[[3]], type = cor_type[3], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "dotplot", paste(cor_type[4], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.mf[[4]], type = cor_type[4], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "dotplot", paste(cor_type[5], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.mf[[5]], type = cor_type[5], annot_type = "MF")
dev.off()

tiff(paste(folder, "annot_bp", "dotplot", paste(cor_type[1], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.bp[[1]], type = cor_type[1], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "dotplot", paste(cor_type[2], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.bp[[2]], type = cor_type[2], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "dotplot", paste(cor_type[3], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.bp[[3]], type = cor_type[3], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "dotplot", paste(cor_type[4], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.bp[[4]], type = cor_type[4], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "dotplot", paste(cor_type[5], ".tiff", sep=''), sep='/'), units="in", width=12, height = 16, res=200)
make_dotplot(all_res.bp[[5]], type = cor_type[5], annot_type = "BP")
dev.off()

#EMAP PLOT
tiff(paste(folder, "annot_mf", "emap", paste(cor_type[1], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.mf[[1]], type = cor_type[1], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "emap", paste(cor_type[2], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.mf[[2]], type = cor_type[2], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "emap", paste(cor_type[3], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.mf[[3]], type = cor_type[3], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "emap", paste(cor_type[4], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.mf[[4]], type = cor_type[4], annot_type = "MF")
dev.off()
tiff(paste(folder, "annot_mf", "emap", paste(cor_type[5], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.mf[[5]], type = cor_type[5], annot_type = "MF")
dev.off()

tiff(paste(folder, "annot_bp", "emap", paste(cor_type[1], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.bp[[1]], type = cor_type[1], annot_type = "BP")
dev.off()
tiff("temp", units="in", width=20, height = 20, res=200)
make_emapplot(all_res.bp[[2]], type = cor_type[2], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "emap", paste(cor_type[3], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.bp[[3]], type = cor_type[3], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "emap", paste(cor_type[4], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.bp[[4]], type = cor_type[4], annot_type = "BP")
dev.off()
tiff(paste(folder, "annot_bp", "emap", paste(cor_type[5], ".tiff", sep=''), sep='/'), units="in", width=20, height = 20, res=200)
make_emapplot(all_res.bp[[5]], type = cor_type[5], annot_type = "BP")
dev.off()








typeeeeeeee

#=================================== CATEGORICAL SCALE ===================================
pearcor_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.pearcor_sig_joined_genes.RDS")
spearcor_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.spearcor_sig_joined_genes.RDS")
pearcor_p_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.pearcor_p_sig_joined_genes.RDS")
spearcor_p_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.spearcor_p_sig_joined_genes.RDS")
randcor_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.randcor_sig_joined_genes.RDS")
# bootcor_fisher_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.bootcor_fisher_sig_joined_genes.RDS")
fisher_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.fisher_sig_joined_genes.RDS")
chisq_res_list = readRDS("all_res_cat3/tissue_spec/all_res_list.chisq_sig_joined_genes.RDS")

pearcor_res.mf = prepare_ck(pearcor_res_list, option = "MF")
spearcor_res.mf = prepare_ck(spearcor_res_list, option = "MF")
pearcor_p_res.mf = prepare_ck(pearcor_p_res_list, option = "MF")
spearcor_p_res.mf = prepare_ck(spearcor_p_res_list, option = "MF")
randcor_res.mf = prepare_ck(randcor_res_list, option = "MF")
# bootcor_fisher_res.mf = prepare_ck(bootcor_fisher_res_list, option = "MF")
fisher_res.mf = prepare_ck(fisher_res_list, option = "MF") #no enrichment
chisq_res.mf = prepare_ck(chisq_res_list, option = "MF") #no enrichment
all_res.mf = list(pearcor_res.mf, spearcor_res.mf, pearcor_p_res.mf, spearcor_p_res.mf,
                  randcor_res.mf, bootcor_fisher_res.mf, fisher_res.mf, chisq_res.mf)
  
pearcor_res.bp = prepare_ck(pearcor_res_list, option = "BP")
spearcor_res.bp = prepare_ck(spearcor_res_list, option = "BP")
pearcor_p_res.bp = prepare_ck(pearcor_p_res_list, option = "BP")
spearcor_p_res.bp = prepare_ck(spearcor_p_res_list, option = "BP")
randcor_res.bp = prepare_ck(randcor_res_list, option = "BP")
# bootcor_fisher_res.bp = prepare_ck(bootcor_fisher_res_list, option = "BP")
fisher_res.bp = prepare_ck(fisher_res_list, option = "BP") 
chisq_res.bp = prepare_ck(chisq_res_list, option = "BP")
all_res.bp = list(pearcor_res.bp, spearcor_res.bp, pearcor_p_res.bp, spearcor_p_res.bp,
                  randcor_res.bp, bootcor_fisher_res.bp, fisher_res.bp, chisq_res.bp)
cor_type = c("pearcor", "spearcor", "pearcor_p", "spearcor_p",
             "randcor", "bootcor_fisher", "fisher", "chisq")

#DOT PLOT
tiff(paste(cor_type[1], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[1]], type = cor_type[1], annot_type = "MF")
dev.off()
tiff(paste(cor_type[2], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[2]], type = cor_type[2], annot_type = "MF")
dev.off()
tiff(paste(cor_type[3], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[3]], type = cor_type[3], annot_type = "MF")
dev.off()
tiff(paste(cor_type[4], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[4]], type = cor_type[4], annot_type = "MF")
dev.off()
tiff(paste(cor_type[5], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[5]], type = cor_type[5], annot_type = "MF")
dev.off()
tiff(paste(cor_type[6], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[6]], type = cor_type[6], annot_type = "MF")
dev.off()
tiff(paste(cor_type[7], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[7]], type = cor_type[7], annot_type = "MF")
dev.off()
tiff(paste(cor_type[8], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.mf[[8]], type = cor_type[8], annot_type = "MF")
dev.off()

tiff(paste(cor_type[1], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[1]], type = cor_type[1], annot_type = "BP")
dev.off()
tiff(paste(cor_type[2], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[2]], type = cor_type[2], annot_type = "BP")
dev.off()
tiff(paste(cor_type[3], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[3]], type = cor_type[3], annot_type = "BP")
dev.off()
tiff(paste(cor_type[4], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[4]], type = cor_type[4], annot_type = "BP")
dev.off()
tiff(paste(cor_type[5], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[5]], type = cor_type[5], annot_type = "BP")
dev.off()
tiff(paste(cor_type[6], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[6]], type = cor_type[6], annot_type = "BP")
dev.off()
tiff(paste(cor_type[7], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[7]], type = cor_type[7], annot_type = "BP")
dev.off()
tiff(paste(cor_type[8], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_dotplot(all_res.bp[[8]], type = cor_type[8], annot_type = "BP")
dev.off()

#EMAP PLOT
tiff(paste(cor_type[1], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[1]], type = cor_type[1], annot_type = "MF")
dev.off()
tiff(paste(cor_type[2], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[2]], type = cor_type[2], annot_type = "MF")
dev.off()
tiff(paste(cor_type[3], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[3]], type = cor_type[3], annot_type = "MF")
dev.off()
tiff(paste(cor_type[4], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[4]], type = cor_type[4], annot_type = "MF")
dev.off()
tiff(paste(cor_type[5], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[5]], type = cor_type[5], annot_type = "MF")
dev.off()
tiff(paste(cor_type[6], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[6]], type = cor_type[6], annot_type = "MF")
dev.off()
tiff(paste(cor_type[7], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[7]], type = cor_type[7], annot_type = "MF")
dev.off()
tiff(paste(cor_type[8], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.mf[[8]], type = cor_type[8], annot_type = "MF")
dev.off()

tiff(paste(cor_type[1], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[1]], type = cor_type[1], annot_type = "BP")
dev.off()
tiff(paste(cor_type[2], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[2]], type = cor_type[2], annot_type = "BP")
dev.off()
tiff(paste(cor_type[3], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[3]], type = cor_type[3], annot_type = "BP")
dev.off()
tiff(paste(cor_type[4], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[4]], type = cor_type[4], annot_type = "BP")
dev.off()
tiff(paste(cor_type[5], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[5]], type = cor_type[5], annot_type = "BP")
dev.off()
tiff(paste(cor_type[6], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[6]], type = cor_type[6], annot_type = "BP")
dev.off()
tiff(paste(cor_type[7], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[7]], type = cor_type[7], annot_type = "BP")
dev.off()
tiff(paste(cor_type[8], ".tiff", sep=''), units="in", width=12, height = 20, res=200)
make_emapplot(all_res.bp[[8]], type = cor_type[8], annot_type = "BP")
dev.off()