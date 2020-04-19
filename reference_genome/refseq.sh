printf "Step 1: Convert transcript ID to gene ID\n"
python GTF_postprocess.py reference.flattened.gtf reference.flattened.ens.gtf reference.genelist.txt
printf "Step 2: Map intron and promoter\n"
python intron_promoter_mapping.py reference.flattened.ens.gtf promoter.gtf intron.gtf
printf "Step 3: Replace „transcripts“ == „transcript_id“\n"
sed 's/transcripts/transcript_id/g' reference.flattened.ens.gtf > temp.txt && mv temp.txt reference.flattened.ens.gtf 
printf "Step 4: Merge intron+promoter+gene+exon together (whole_genome.gtf)\n"
cat promoter.gtf intron.gtf reference.flattened.ens.gtf > whole_genome.gtf && bedtools sort -i whole_genome.gtf > whole_genome_sorted.gtf && mv whole_genome_sorted.gtf whole_genome.gtf
printf "Step 5: Get TSS\n" 
awk '{if (($0 ~ "aggregate_gene") && ($7=="+")) print $1"\t"$2"\tTSS\t"$4"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10; else if (($0 ~ "aggregate_gene") && ($7=="-")) print $1"\t"$2"\tTSS\t"$5"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10}' reference.flattened.ens.gtf > TSS.gtf
printf "Step 6: Merge intron+promoter+gene+exon+TSS (final_wholegenome.gtf)\n"
cat whole_genome.gtf TSS.gtf > temp.gtf && bedtools sort -i temp.gtf > final_wholegenome.gtf && rm temp.gtf