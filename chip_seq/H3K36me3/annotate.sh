cat /dev/null > annotate_manorm.log
echo "Start Time: $(date)" >> annotate_manorm.log
mkdir normalizedcounts
mkdir annotatedcounts

echo "===> Begin annotating"
#Get chr, start, end, M-value, p-value, normedcount1, normedcount2
echo "===> 1: Get chr, start, end, M-value, p-value, normedcount1, normedcount2"
for x in manorm_result/E*/*.xls
do (
    echo $x
    f=${x##*/}
    awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$7"\t"$9"\t"$10}' $x > normalizedcounts/$f.bed
) &
done
wait

#annotate all xls.bed in normalizedcounts into annotatedcounts
echo "===> 2: Annotate all xls.bed in normalizedcounts into annotatedcounts"
annotate_manorm_parallel_exon() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[5]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[6]}')
    if (ls normalizedcounts| grep "$epi1" ) && (ls normalizedcounts | grep "$epi2" )
    then (
        if (ls annotatedcounts | grep "exon_"$epi1"_"$epi2"")
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_manorm.log
        )
        else (
            echo "Annotating pair "$epi1"_"$epi2" - EXON" >> annotate_manorm.log
            REF_GEN_EXON=$REFGEN/reference_genome.gtf
            HIS_COUNT=normalizedcounts/"$epi1"_"$epi2"_all_MAvalues.xls.bed
            ANNOT_HIS_COUNT=annotatedcounts/exon_"$epi1"_"$epi2".txt
            bedtools intersect -a $REF_GEN_EXON -b $HIS_COUNT -wo -loj -bed > $ANNOT_HIS_COUNT
            )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_manorm.log
    fi
}
annotate_manorm_parallel_pi() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[5]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[6]}')
    if (ls normalizedcounts| grep "$epi1" ) && (ls normalizedcounts | grep "$epi2" )
    then (
        if (ls annotatedcounts | grep "pi_"$epi1"_"$epi2"")
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_manorm.log
        )
        else (
            echo "Annotating pair "$epi1"_"$epi2" - PI" >> annotate_manorm.log
            REF_GEN_PI=$REFGEN/reference_genome.promoter_intron.gtf
            HIS_COUNT=normalizedcounts/"$epi1"_"$epi2"_all_MAvalues.xls.bed
            ANNOT_HIS_COUNT=annotatedcounts/pi_"$epi1"_"$epi2".txt
            bedtools intersect -a $REF_GEN_PI -b $HIS_COUNT -wo -loj -bed > $ANNOT_HIS_COUNT
            )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_manorm.log
    fi
}

for chunk in pair_chunks/*
do
    while read p;
    do
    annotate_manorm_parallel_exon &
    done < $chunk
    wait
    echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - EXON" >> annotate_manorm.log
done
echo "End Time: $(date)" >> annotate_manorm.log

for chunk in pair_chunks/*
do
    while read p;
    do
    annotate_manorm_parallel_pi &
    done < $chunk
    wait
    echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - PI" >> annotate_manorm.log
done
echo "End Time: $(date)" >> annotate_manorm.log