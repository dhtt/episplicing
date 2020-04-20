echo "===> Begin annotating"
  
#annotate all xls.bed in normalizedcounts into annotatedcounts
echo "===> 1: Annotate all normedratio.csv in normalizedcounts into annotatedcounts"
for f in normalizedcounts/*normedratio.csv
do
    (
    echo $f
    FILENAME="${f#*/}.txt"
    bedtools intersect -a $REFGEN/reference_genome.gtf -b $f -wo -loj -bed > annotatedcounts/$FILENAME
    bedtools intersect -a $REFGEN/reference_genome.promoter_intron.gtf -b $f -wo -loj -bed > annotatedcounts/"pi_$FILENAME"
    grep -v "aggregate_gene" annotatedcounts/$FILENAME > annotatedcounts/"exon_$FILENAME"
    rm annotatedcounts/$FILENAME
    )&
done
wait

#annotate all methdiff25 in diff into annotateddiff
echo "===> 2: Annotate all methdiff25 in diff into annotateddiff"
for f in normalizedcounts/*diff.txt
do
    (
    echo $f
    FILENAME="${f#*/}.txt"
    bedtools intersect -a $REFGEN/reference_genome.gtf -b $f -wo -loj -bed > annotateddiff/$FILENAME
    bedtools intersect -a $REFGEN/reference_genome.promoter_intron.gtf -b $f -wo -loj -bed > annotateddiff/"pi_$FILENAME"
    grep -v "aggregate_gene" annotateddiff/$FILENAME > annotateddiff/"exon_$FILENAME"
    rm annotateddiff/$FILENAME
    )&

done
wait
echo "===> Finished annotating"
