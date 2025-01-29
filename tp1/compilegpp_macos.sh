#!/usr/bin/env bash

echo -e "\033[33mWarning: You need to install XQuartz.\033[0m"

file="$1"

# Make sure the file exists
if [ ! -f "$file" ]; then
    echo "File does not exist"
    exit 1
fi

# Make sure the file is a cpp file
if [[ "$file" != *.c ]]; then
    echo "File is not a cpp file"
    exit 1
fi

file_no_ext="$(basename $file .c)"
echo "Compiling ${file}"

# Compile
g++ -o "$file_no_ext" "$file" -I/usr/X11/include -L/usr/X11/lib -lm -lX11

# Get output code of the compilation
output=$?
if [ $output -eq 0 ]; then
    echo "Saved executable at ./${file_no_ext}"
else
    echo "Compilation failed"
    exit 1
fi
