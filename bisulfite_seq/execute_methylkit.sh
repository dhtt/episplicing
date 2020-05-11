cat /dev/null > execute_methylkit.log
echo "Start Time: $(date)" >> execute_methylkit.log
echo "===> START METHYLKIT-TING ALL PAIRS"
f="$MET"summary

execute_methylkit_parallel() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls $f | grep "$epi1" ) && (ls $f | grep "$epi2" )
    then (
        if (ls "$MET"normalizedcounts | grep "$epi1"_"$epi2" )
        then (
            echo "Pair "$epi1"_"$epi2" already exists" >> execute_methylkit.log
        )
        else (
            Rscript "$MET"methylKit_analysis.R -a "$epi1" -b "$epi2" -f  "$f" -c 10 -n 8
            echo "Methylkit-ing pair "$epi1"_"$epi2"" >> execute_methylkit.log
        )
        fi
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> execute_methylkit.log
    fi
}

for chunk in pair_chunks/*
do
    while read p;
    do
    execute_methylkit_parallel &
    done < $chunk
    wait
    echo "===> FINISHED METHYLKIT-ING ALL PAIRS IN "$chunk"" >> execute_methylkit.log
done
echo "End Time: $(date)" >> execute_methylkit.log
