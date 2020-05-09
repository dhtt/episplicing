cat /dev/null > execute_dexseq.log
echo "===> START DEXSEQ-ING ALL PAIRS"

execute_dexseq_parallel() {
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls . | grep "$epi1" ) && (ls . | grep "$epi2" )
    then (
        Rscript "$EXP"DEXSeq_analysis.R -a "$epi1" -b "$epi2" -f  "$EXP"backupcount_NCBI/ -g $REFGEN/reference_genome.gtf -n 5
        echo "DEXseq-ing pair "$epi1"_"$epi2"" >> execute_dexseq.log
        )
    else echo "Pair "$epi1"_"$epi2" does not exist" >> execute_dexseq.log
    fi
}

while read p; 
do
    execute_dexseq_parallel &
done < "$EXP"all_pairs.txt
wait
echo "===> FINISHED DEXSEQ-ING ALL PAIRS"

