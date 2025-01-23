#!/usr/bin/env bash

file="$1"

# Make sure the file exists
if [ ! -f "$file" ]; then
    echo "File does not exist"
    exit 1
fi

# Make sure the file is a cpp file
if [[ "$file" != *.cpp ]]; then
    echo "File is not a cpp file"
    exit 1
fi

file_no_ext="$(basename $1 .cpp)"
echo "Compiling ${file}"

# Get output code of the following command
g++ -o "$file_no_ext" "$file" -lm -lX11

output=$?
if [ $output -eq 0 ]; then
    echo "Saved executable at ./${file_no_ext}"
else
    echo "Compilation failed"
    exit 1
fi
