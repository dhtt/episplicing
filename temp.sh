mkdir subset
cp normalizedcounts/E003_E004_all_MAvalues.xls.bed normalizedcounts/E003_E100_all_MAvalues.xls.bed  normalizedcounts/E004_E100_all_MAvalues.xls.bed normalizedcounts/E094_E095_all_MAvalues.xls.bed normalizedcounts/E003_E005_all_MAvalues.xls.bed normalizedcounts/E004_E005_all_MAvalues.xls.bed normalizedcounts/E005_E094_all_MAvalues.xls.bed normalizedcounts/E094_E100_all_MAvalues.xls.bed normalizedcounts/E003_E094_all_MAvalues.xls.bed normalizedcounts/E004_E094_all_MAvalues.xls.bed normalizedcounts/E005_E095_all_MAvalues.xls.bed normalizedcounts/E095_E100_all_MAvalues.xls.bed normalizedcounts/E003_E095_all_MAvalues.xls.bed normalizedcounts/E004_E095_all_MAvalues.xls.bed normalizedcounts/E005_E100_all_MAvalues.xls.bed subset/
cd subset/

for f in *MAvalues.xls.bed
do (
    echo $f
    bedtools intersect -a $REFGEN/reference_genome.fl200.gtf -b $f -wo -loj -bed > fl200_$f.txt 
) &
done
wait 

for f in fl200_*
do (
    echo $f
    bedtools groupby -i $f -g 1-9 -c 13,14 -o mode > fl200_$f
) &
done
wait
echo "Finish"

mkdir subset
mv fl200_fl200* subset/
tar -cvf subset.tar.gz subset
