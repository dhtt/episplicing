echo "===> Begin annotating"
for d in E*/ #get all xls out in normalizedcounts folder
do
    echo "$d"
    for x in $d/*.xls; do awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$10}' $x > $x.bed; done
    cp $d/*.xls.bed normalizedcounts
done

for f in normalizedcounts/*.bed #annotate all xls.bed in normalizedcounts into annotatedcounts
do  
    (echo $f
    FILENAME="${f#*/}.txt"
    bedtools intersect -a $REFGEN/reference_genome.gtf -b $f -wo -loj -bed > annotatedcounts/$FILENAME 
    grep -v "aggregate_gene" annotatedcounts/$FILENAME > annotatedcounts/"exon_$FILENAME") &
    
done 
wait
echo "===> Finished annotating"