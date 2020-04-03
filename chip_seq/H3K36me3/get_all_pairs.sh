ls *.narrowPeak *tagAlign > filenames.txt
cat filenames.txt | paste -sd '\t\n' > filelist.txt
#rm filenames.txt

awk '{ a[NR]=$0 }
       END{ for(i=1;i<=NR;i++)
              for(j=i+1;j<=NR;j++)
                print a[i], a[j] }' filelist.txt > all_pairs.txt

