echo "=====> Retrieving histone files"
# TASK1 =====================================
echo ".....> Retrieving links"
lynx -dump https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/ | awk '/http/{print $2}' > narrowpeak_links.txt
lynx -dump https://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated/ | awk '/http/{print $2}' > chip_links.txt

# TASK2 =====================================
echo ".....> Downloading from links"
types=("H3K4me1" "H3K4me3" "H3K9me3" "H3K9ac" "H3K27me3" "H3K27ac" "H3K36me3")
for i in ${types[*]}
do
    type="$i"
    FILE_PATH="/home/dhthutrang/files/chip_seq/$type"
    FILE_NAME="$type"_links.txt

    grep $type narrowpeak_links.txt >> $FILE_NAME
    grep $type chip_links.txt >> $FILE_NAME

    echo $FILE_PATH > temp1.txt
    grep -hi -f epigenomes.txt $FILE_NAME >> temp1.txt
    grep -v ".tbi" temp1.txt > temp2.txt

    mv temp2.txt $FILE_NAME
    rm temp*
done

for file in H*links.txt 
do
    ./aria2files.sh $file
done

# TASK3 =====================================
echo ".....> Unzipping narrowpeak.gz and tagalign.gz"
unzip_gz(){
    type="$i"
    FILE_PATH="/home/dhthutrang/files/chip_seq/$type"
    echo $FILE_PATH
    cd $FILE_PATH && gunzip *
}
types=("H3K4me1" "H3K4me3" "H3K9me3" "H3K9ac" "H3K27me3" "H3K27ac" "H3K36me3")
for i in ${types[*]}
do
    unzip_gz $i &
done

# TASK3 SUB =====================================
echo ".....> Treat gene name incompatibility in DEXSeq"
for i in *count.txt 
do 
    sed -i 's/Em:AC005003.4/EmAC005003.4/g' $i
    sed -i 's/Em:AC008101.5/EmAC008101.5/g' $i
done

