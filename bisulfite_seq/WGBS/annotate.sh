cat /dev/null > annotate_methylkit.log
echo "Start Time: $(date)" >> annotate_methylkit.log

annotate_methylkit_parallel_exon() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls summary | grep "$epi1" ) && (ls summary | grep "$epi2" )
    then (
        if (ls annotatedcounts | grep "exon_"$epi1"_"$epi2"" )
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_methylkit.log
        )
        else (
            echo "Annotating pair "$epi1"_"$epi2" - EXON" >> annotate_methylkit.log
            REF_GEN_EXON=$REFGEN/reference_genome.gtf
            MET_COUNT=summary/"$epi1"_"$epi2"_normedratio.csv
            ANNOT_MET_COUNT=annotatedcounts/exon_"$epi1"_"$epi2"_normedratio.csv.txt
            bedtools intersect -a $REF_GEN_EXON -b $MET_COUNT -wo -loj -bed > $ANNOT_MET_COUNT
             )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_methylkit.log
    fi
}
annotate_methylkit_parallel_pi() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls summary | grep "$epi1" ) && (ls summary | grep "$epi2" )
    then (
        if (ls annotatedcounts | grep "pi_"$epi1"_"$epi2"" )
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_methylkit.log
        )
        else (
            echo "Annotating pair "$epi1"_"$epi2" - PI" >> annotate_methylkit.log
            REF_GEN_PI=$REFGEN/reference_genome.promoter_intron.gtf
            MET_COUNT=summary/"$epi1"_"$epi2"_normedratio.csv
            ANNOT_MET_COUNT=annotatedcounts/pi_"$epi1"_"$epi2"_normedratio.csv.txt
            bedtools intersect -a $REF_GEN_PI -b $MET_COUNT -wo -loj -bed > $ANNOT_MET_COUNT
        )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_methylkit.log
    fi
}

for chunk in pair_chunks/*
do
    while read p;
    do
    annotate_methylkit_parallel_exon &
    done < $chunk
    wait
    echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - EXON" >> annotate_methylkit.log
done
echo "End Time: $(date)" >> annotate_methylkit.log

for chunk in pair_chunks/*
do
    while read p;
    do
    annotate_methylkit_parallel_pi &
    done < $chunk
    wait
    echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - PI" >> annotate_methylkit.log
done
echo "End Time: $(date)" >> annotate_methylkit.log

cat /dev/null > annotate_methylkit_diff.log
echo "Start Time: $(date)" >> annotate_methylkit_diff.log

annotate_methylkit_diff_parallel() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls summary | grep "$epi1" ) && (ls summary | grep "$epi2" )
    then (
        if (ls annotateddiff | grep "exon_"$epi1"_"$epi2"" )
        then (
            echo "Pair "$epi1"_"$epi2" already annotated" >> annotate_methylkit_diff.log
        )
        else (
            echo "Annotating pair "$epi1"_"$epi2" - EXON" >> annotate_methylkit_diff.log
            REF_GEN=$REFGEN/reference_genome.whole.sorted.gtf
            MET_DIFF=summary/"$epi1"_"$epi2"_diff.txt
            ANNOT_MET_DIFF=annotateddiff/"$epi1"_"$epi2"_diff.txt.txt
            bedtools intersect -a $REF_GEN -b $MET_DIFF -wo -loj -bed > $ANNOT_MET_DIFF
            )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> annotate_methylkit_diff.log
    fi
}
for chunk in pair_chunks/*
do
    while read p;
    do
    annotate_methylkit_diff_parallel &
    done < $chunk
    wait
    echo "===> FINISHED ANNOTATING ALL PAIRS IN "$chunk" - DIFF" >> annotate_methylkit_diff.log
done
echo "End Time: $(date)" >> annotate_methylkit_diff.log
