reference_genome = data.frame(fread("/Users/dhthutrang/Documents/BIOINFO/Episplicing/episplicing/mrna_seq/reference_genome.gtf", header=FALSE))
reference_genome = as.data.table(reference_genome[grep("exonic_part", reference_genome$V3),9])
temp = as.data.table(str_split_fixed(reference_genome$V1, '; ', 3))
temp1 = as.data.frame(gsub('transcript_id ', '', temp$V2))
temp2 = as.data.frame(gsub('gene_id ', '', temp$V1))
temp$V2 = temp1
temp$V1 = temp2

temp3 <- temp %>%
  group_by(V1) %>%
  mutate(transcript = paste(V2, collapse = '+'))

temp3 <- temp3 %>%
  group_by(V1,transcript) %>%
  select(V1, transcript) %>%
  unique()
temp4 = strsplit(temp3$transcript, "+", fixed=TRUE)
temp4 = lapply(temp4, function(x) gsub('\"', "", x))
temp5 = lapply(temp4, unique)
temp5 = lapply(temp5, function(x) paste(x, collapse = ','))
temp6 = transpose(as.data.table(temp5))
all_transcripts = cbind(as.data.frame(temp3$V1), temp6$V1)
colnames(all_transcripts) = c("gene", "transcript")
all_transcripts$gene = gsub("\"", "", all_transcripts$gene)
all_transcripts <- all_transcripts %>%
  mutate(n_transcript = str_count(transcript, ',') + 1 ) %>%
  mutate(n_gene = str_count(gene, '\\+') + 1 )
constitutive_genes <- all_transcripts %>%
  filter(n_transcript == 1)
as_genes <- all_transcripts %>%
  filter(n_transcript > 1)

fwrite(constitutive_genes, file= "constitutive_genes.txt",quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
fwrite(as_genes, file= "as_genes.txt",quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)










