#!/bin/bash

# Check if both files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <file1> <file2>"
    exit 1
fi

file1="$1"
file2="$2"

# Check if both files exist
if [ ! -f "$file1" ]; then
    echo "Error: File '$file1' not found!"
    exit 1
fi

if [ ! -f "$file2" ]; then
    echo "Error: File '$file2' not found!"
    exit 1
fi

# Compare the files using diff
echo "Comparing $file1 and $file2..."
diff_output=$(diff -u "$file1" "$file2")

if [ -z "$diff_output" ]; then
    echo "The files are identical."
else
    echo "Differences found:"
    echo "$diff_output"
fi
