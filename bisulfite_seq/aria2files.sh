#!/bin/bash
#author: hmj6jmh/https://askubuntu.com/questions/214018/how-to-make-wget-faster-or-multithreading

filename="$1" # get filename from command line argument

while read -r line
do
    if [ "$line" ] # skip blank lines
    then
        if [[ "$line" =~ (https?|ftp)\:\/\/ ]] # line contains a URL, download file
        then
            echo "URL: '$line'"
            aria2c --file-allocation=none -c -x 10 -s 10 -d "$currdir" "$line"
        else # line contains a directory name, create directory if not already present
            echo "Directory: '$line'"
            currdir="$line"
            if [ ! -d "$currdir" ]
            then
                mkdir -p "$currdir" # '-p' enables creation of nested directories in one command
            fi
        fi
    fi
done < "$filename"
