echo "=====> Begin mRNA preprocessing pipeline"
cd "$(dirname "$0")"

if (test -f mRNA_seq_BEDaccessionnumber.txt) && (test -f hg19.genome)
then
    echo "...File requirement met."
else
    echo "Error: check mRNA_seq_BEDaccessionnumber.txt or hg19.genome"
    exit 1
fi

# TASK1 =====================================
unzip_bed(){
    ACC_NO=${f%%_*}
    echo "unzip_bed $ACC_NO" >> log_code1.txt
    if ! grep -Fxq "$ACC_NO" mRNA_seq_BEDaccessionnumber.txt
    then 
        echo "Err: Cannot find $ACC_NO" >> log_code1.txt
    else
        gunzip $f || echo "Err: Skipping file $f" >> log_code1.txt
    fi
}
echo ".....> Unzipping bed.gz"
for f in *.bed.gz
do
    unzip_bed $f &
done
wait

# TASK1 - SUB1 =====================================
echo ".....> Check if BED has 6 columns"
for f in *.bed
do
    ACC_NO=${f%%_*}
    cat $f | awk -v acc_no=$ACC_NO '{if (NF!=6) print acc_no; exit}' >> check.txt
done

check_col(){
    ACC_NO=${f%%_*}
    if grep -Fxq "$ACC_NO" check.txt
    then
        echo $ACC_NO
        awk 'BEGIN {OFS="\t"}{$6=$5; $5=0; print $0}' $f > "$ACC_NO"_extracol.bed
        rm $f
    fi
}
for f in *.bed
do
    check_col $f &
done

# TASK2 =====================================
bed_to_bam(){
    ACC_NO=${f%%_*}
    echo "bed_to_bam $ACC_NO" >> log_code1.txt
    (bedtools bedtobam -i $f -g hg19.genome > $ACC_NO.bam && rm $f) || echo "Err: Cannot convert BED to BAM $ACC_NO" >> log_code1.txt
}
echo ".....> Converting BED to BAM"
for f in *.bed
do 
    bed_to_bam $f &
done
wait


# TASK3 =====================================
bam_to_sam(){
    ACC_NO=${f%%.*}
    echo "bam_to_sam $ACC_NO" >> log_code1.txt
    (samtools sort -n -O SAM $f > $ACC_NO.sam) || echo "Err: Cannot convert BAM to SAM $ACC_NO" >> log_code1.txt
}
echo ".....> Converting BAM to SAM"
for f in *.bam
do 
    bam_to_sam $f &
done
wait

echo "=====> Finished mRNA pipeline"