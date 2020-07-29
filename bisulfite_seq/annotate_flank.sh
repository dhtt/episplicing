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
            REF_GEN=$REFGEN/reference_genome.fl200.gtf
            MET_DIFF=summary/"$epi1"_"$epi2"_diff.txt
            ANNOT_MET_DIFF=$FLANK/met/"$epi1"_"$epi2"_diff.txt.txt
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


for f in $FLANK/met/*
do (
    echo $f 
    bedtools groupby -i $f -g 1-9 -c 13,14 -o mode > $f.fl.txt
) &
done
wait
echo "Finish"

cd $FLANK/met
rm ./!(*txt.fl.txt*)