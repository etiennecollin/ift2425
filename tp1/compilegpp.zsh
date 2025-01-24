#!/usr/bin/env zsh

local usage_message="Usage:\n\tTo run: ./compilegpp.sh -t|--target <target> -f|--file <file.cpp>\n\tTo list targets: ./compilegpp.sh -l|--list"
local target_list="linux\nmacos"
local target=""
local file=""
local gpp_arguments=""
local arguments_pair="$(($# % 2))"

# Check if no arguments are provided
if [ $# -eq 0 ]; then
    echo "$usage_message"
    exit 1
fi

# Make sure the arguments are even or there is only one argument

# Parse the arguments
local i=0
while [ $# -gt 0 ]; do
    case $1 in
    -t | --target)
        if [ "$arguments_pair" -ne 0 ]; then
            echo "Error: Missing value for a flag."
            echo "$usage_message"
            exit 1
        fi
        target="$2"
        shift 2
        ;;
    -f | --file)
        if [ "$arguments_pair" -ne 0 ]; then
            echo "Error: Missing value for a flag."
            echo "$usage_message"
            exit 1
        fi
        file="$2"
        shift 2
        ;;
    -l | --list)
        echo "$target_list"
        exit 0
        ;;
    *)
        echo "$usage_message"
        exit 1
        ;;
    esac
    i=$((i + 1))
done

# Check if target is linux or macos
if [ "$target" = "linux" ]; then
    gpp_arguments="-lm -lX11"
elif [ "$target" = "macos" ]; then
    gpp_arguments="-I/usr/X11/include -L/usr/X11/lib -lm -lX11"
    echo "Warning: You need to install XQuartz and run the compilation in xterm."
else
    echo "Error: Invalid target."
    echo "$usage_message"
    exit 1
fi
# Format the string as a list of args
gpp_arguments=( ${(z)gpp_arguments} )

# Make sure the file is provided and exists
if [ "$file" = "" ] || [ ! -f "$file" ]; then
    echo "Error: File does not exist"
    exit 1
fi

# Make sure the file is a cpp file
if [[ "$file" != *.c ]]; then
    echo "Error: File is not a cpp file"
    exit 1
fi

file_no_ext="$(basename $file .c)"
echo "Compiling ${file} for ${target}..."

# Get output code of the following command
g++ -o "$file_no_ext" "$file" $gpp_arguments

output=$?
if [ $output -eq 0 ]; then
    echo "Saved executable at ./${file_no_ext}"
else
    echo "Compilation failed"
    exit 1
fi
