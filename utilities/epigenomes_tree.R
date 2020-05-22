#===== LOAD PACKAGES ======
library("data.table", quietly=TRUE)
library("dplyr", quietly=TRUE)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 

d1 = fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/utilities/epigenomes.txt",
           header = TRUE, sep = '\t')
vertices = data.frame(name = unique(c(as.character(d1$from), as.character(d1$to)))) 
vertices$group = d1$from[ match( vertices$name, d1$to) ]
vertices$id=NA
myleaves=which(is.na( match(vertices$name, d1$from) ))
nleaves=length(myleaves)
vertices$id[ myleaves ] = seq(1:nleaves)
mygraph <- graph_from_data_frame( d1, vertices=vertices )

tiff("epigenomes.tiff", units="in", width=7.6, height=6.2, res=300)
ggraph(mygraph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal(colour="grey") +
  geom_node_text(aes(x=x-0.3, y=y, label=name, filter=leaf==FALSE, hjust = 1, colour=group), 
                 size=2.8, alpha=1) +
  geom_node_text(aes(x=x, y=y-0.15, label=name, filter=leaf, angle = 30, hjust = 1, colour=group), 
                 size=2.8, alpha=1) +
  geom_node_point(aes(x = x, y=y, colour=group, alpha=0.2)) +
  scale_colour_manual(values = brewer.pal(8, "Dark2")) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"inches"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
dev.off()








