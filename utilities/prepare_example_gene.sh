if [[ $# -ne 3 ]]
then 
    echo "Please input epi_id1, epi_id2 and gene name"
    exit 1
fi
mkdir example_gene
epi1=$1
epi2=$2
gene=$3

#Prepare exp file
echo "====> Prepare exp file"
grep "\"$gene\"" $REFGEN/hg19.ncbiRefSeq.gtf > example_gene/"$gene"_NCBI.gtf
sed -i 's/_id//g' example_gene/"$gene"_NCBI.gtf

grep "groupID" $EXP/backupcount_NCBI/res/"$epi1"_"$epi2"_res.csv > temp1.txt
grep $gene$'\t' $EXP/backupcount_NCBI/res/"$epi1"_"$epi2"_res.csv >> temp1.txt
grep "\"$gene\"" $REFGEN/reference_genome.gtf > temp2.txt
paste -d '\t' temp1.txt temp2.txt > example_gene/"$epi1"_"$epi2"_res.csv
rm temp1.txt temp2.txt

#Prepare met file
echo "====> Prepare met file"
##Met diff
grep "\"$gene\"" $MET/annotateddiff/"$epi1"_"$epi2"_diff.txt.txt > example_gene/"$epi1"_"$epi2"_diff.txt.txt
##Met count
grep "\"$gene\"" $MET/annotatedcounts/exon_"$epi1"_"$epi2"_normedratio.csv.txt > example_gene/exon_"$epi1"_"$epi2"_normedratio.csv.txt
grep "\"$gene\"" $MET/annotatedcounts/pi_"$epi1"_"$epi2"_normedratio.csv.txt > example_gene/pi_"$epi1"_"$epi2"_normedratio.csv.txt

#Prepare his file
echo "====> Prepare his file"
for his in $HIS/H*/
do (
    temp=${his%*/}
    type=${temp##*/}
    grep "\"$gene\"" $HIS/$type/annotatedcounts/exon_"$epi1"_"$epi2".txt > $UTI/example_gene/exon_"$epi1"_"$epi2"."$type".txt
    grep "\"$gene\"" $HIS/$type/annotatedcounts/pi_"$epi1"_"$epi2".txt > $UTI/example_gene/pi_"$epi1"_"$epi2"."$type".txt
)
done
#Finishing
tar -cvf example_"$epi1"_"$epi2"_"$gene".tar.gz example_gene
rm -r example_gene

epi1=$1
epi2=$2
gene=$3
grep "\"$gene\"\|"groupID"" "$epi1"_"$epi2"_res.csv 
