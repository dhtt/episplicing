ls *.narrowPeak *tagAlign > filenames.txt
cat filenames.txt | paste -sd '\t\n' > filelist.txt
rm filenames.txt

awk '{ a[NR]=$0 }
       END{ for(i=1;i<=NR;i++)
              for(j=i+1;j<=NR;j++)
                print a[i], a[j] }' filelist.txt > all_pairs.temp1.txt

cp $EXP/all_pairs.txt ./all_pairs.temp2.txt

while read p 
do (
    epi1=$(echo $p | awk '{split($0, a, " "); print a[1]}')
    epi2=$(echo $p | awk '{split($0, a, " "); print a[2]}')
    if (ls . | grep "$epi1" ) && (ls $f | grep "$epi2" )
    then
    (
    files=$(cat all_pairs.temp1.txt | grep "$epi1" | grep "$epi2")
    echo $files"        "$epi1" "$epi2"" >> all_pairs.txt
    )
    fi
)
done < all_pairs.temp2.txt

mkdir pair_chunks
split -l 11 --numeric-suffixes all_pairs.txt pair_chunks/chunk_

rm all_pairs.temp1.txt all_pairs.temp2.txt