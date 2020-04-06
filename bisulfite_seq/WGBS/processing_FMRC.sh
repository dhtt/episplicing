#Processing FM and RC from Roadmaps Epigenomes#
##Prepare chromosomes and postions##
awk '{split (FILENAME,name,"."); split ($0,a); print name[1]"\t"a[1]"\t+"}' *.fm > all_chroms.txt

##Separate column in FM and RC into files for selected epigenomes##
while read line
do
    (idx=$(echo $line | awk '{print $1}')
    NAME=$(echo $line | awk '{print $2;}')
    echo $idx
    echo $NAME
    awk -v i="$idx" '{print $i}' FM/*.fm > FM/"$NAME.txt"
    awk -v i="$idx" '{print $i}' RC/*.rc > RC/"$NAME.txt") &
done < epi_ids.txt
wait
echo "===> Finished separating files"

##Generating counts from FM and RC for each epigenome##
for f in FM/E*.txt
do
    NAME="${f#*/}"
    paste -d"\t\t" all_chroms.txt FM/$NAME RC/$NAME > "normalizedcounts"/"$NAME"_"count.txt"
done
wait
echo "===> Finished generating counts for epigenomes"
