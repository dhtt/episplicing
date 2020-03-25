echo "=====> Begin counting"
cd "$(dirname "$0")"

if (test -f mRNA_seq_BEDaccessionnumber.txt) && (test -f reference.flattened.ens.gtf)
then
    echo "...File requirement met."
else
    echo "Err: Check mRNA_seq_BEDaccessionnumber.txt or reference.flattened.ens.gtf"
    exit 1
fi

# COUNTING =====================================
dexseq_count() {
    f="$p".sam
    ACC_NO=${f%%.*}
    echo "counting $ACC_NO" >> log_dexseqcount.txt
    (python ~/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_count.py reference.flattened.ens.gtf $f "$ACC_NO"_count.txt) || (echo "Err: Counting $ACC_NO failed" >> log_dexseqcount.txt)
}

while read p; 
do
    dexseq_count &
done < mRNA_seq_BEDaccessionnumber.txt #temp.txt is accession_number (with order sometimes)

# Change name to contain epigenome ID =====================================
for f in *count.txt
do
    ACC_NO=${f%%_*}
    echo "Old name: $ACC_NO"
    ID=$(grep -ih $ACC_NO temp.txt)
    prefix_ID=${ID%%_*}
    NEWNAME="$prefix_ID"_"$f"
    mv $f $NEWNAME
done