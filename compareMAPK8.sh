#!/bin/bash

# Check the number of arguments
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <file1> <file2> [columns_to_ignore]"
    echo "Example: $0 file1.csv file2.csv 1 3"
    exit 1
fi

# Arguments
file1="$1"
file2="$2"
shift 2  # Shift arguments so that additional arguments become $1, $2, ...

# Check if the files exist
if [[ ! -f "$file1" || ! -f "$file2" ]]; then
    echo "One or both files do not exist. Check the file names and try again."
    exit 1
fi

# If columns to ignore are specified, store them
cut_columns=("$@")  # Remaining arguments are column numbers to be ignored

# Function to filter and sort columns within each row
normalize_file() {
    local file="$1"
    if [[ ${#cut_columns[@]} -eq 0 ]]; then
        awk -F',' '{ 
            split($0, cols, ","); 
            asort(cols); 
            print join(cols, ",") 
        }' "$file" | sort -u
    else
        awk -F',' -v cols="${cut_columns[*]}" '
        BEGIN {split(cols, ignore_cols, " ")}
        {
            out = ""
            for (i = 1; i <= NF; i++) {
                found = 0
                for (j in ignore_cols) {
                    if (i == ignore_cols[j]) {
                        found = 1
                        break
                    }
                }
                if (!found) {
                    out = (out == "" ? $i : out "," $i)
                }
            }
            split(out, fields, ",")
            asort(fields)
            new_out = fields[1]
            for (k = 2; k in fields; k++) {
                new_out = new_out "," fields[k]
            }
            print new_out
        }' "$file" | sort -u
    fi
}

# Filtered and normalized files (excluding specific columns and sorting rows and columns)
filtered_file1=$(mktemp)
filtered_file2=$(mktemp)

normalize_file "$file1" > "$filtered_file1"
normalize_file "$file2" > "$filtered_file2"

# Check if all lines from the first file are present in the second file
missing_lines=$(comm -23 "$filtered_file1" "$filtered_file2")

# Delete temporary files
rm "$filtered_file1" "$filtered_file2"

if [[ -z "$missing_lines" ]]; then
    echo "All lines from $file1 are present in $file2 (ignoring column order and excluding columns: ${cut_columns[*]})."
else
    echo "The following lines from $file1 are missing in $file2 (ignoring column order and excluding columns: ${cut_columns[*]}):"
    echo "$missing_lines"
fi
