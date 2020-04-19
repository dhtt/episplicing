echo "===> Begin annotating"
#Get chr, start, end, M-value, p-value, normedcount1, normedcount2
echo "===> 1: Get chr, start, end, M-value, p-value, normedcount1, normedcount2"
for x in E*/*.xls
do
        awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$7"\t"$9"\t"$10}' $x > $x.bed
done

#get all xls.bed out in normalizedcounts folder
echo "===> 2: Move all xls.bed out in normalizedcounts folder"
for d in E*/
do
    #echo "$d"
    cp $d/*.xls.bed normalizedcounts
    rm $d/*.xls.bed
done

#annotate all xls.bed in normalizedcounts into annotatedcounts
echo "===> 3: Annotate all xls.bed in normalizedcounts into annotatedcounts"
for f in normalizedcounts/*.bed
do  
    (
    #echo $f
    FILENAME="${f#*/}.txt"
    bedtools intersect -a $REFGEN/reference_genome.gtf -b $f -wo -loj -bed > annotatedcounts/$FILENAME 
    bedtools intersect -a $REFGEN/reference_genome.promoter_intron.gtf -b $f -wo -loj -bed > annotatedcounts/"pi_$FILENAME"
    grep -v "aggregate_gene" annotatedcounts/$FILENAME > annotatedcounts/"exon_$FILENAME"
    rm annotatedcounts/$FILENAME
    ) &
    
done 
wait
echo "===> Finished annotating"
