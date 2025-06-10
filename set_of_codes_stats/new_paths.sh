#!/bin/bash

GENES=(
 # "ENSG00000058404"
 # "ENSG00000197122"
 # "ENSG00000140416"
 # "ENSG00000139220"
 # "ENSG00000107643"
 # "ENSG00000187122"
 # "ENSG00000128641"
 # "ENSG00000148737"
 # "ENSG00000010810"
 # "ENSG00000169855"
 # "ENSG00000185008"
 # "ENSG00000092421"
 # "ENSG00000050748"



 "ENSMUSG00000019843"
 "ENSMUSG00000024985"
 "ENSMUSG00000018417"
 "ENSMUSG00000019647"
 "ENSMUSG00000020366"
 "ENSMUSG00000021936"
 "ENSMUSG00000022883"
 "ENSMUSG00000025020"
 "ENSMUSG00000027646"
 "ENSMUSG00000032366"
 "ENSMUSG00000052516"
 "ENSMUSG00000053825"
 "ENSMUSG00000057897"
)


OUTPUT="summary_counters_mousee.txt"
echo -e "Gene\tCount1(%)\tCount2(%)\tCount3(%)\tCount4(%)" > "$OUTPUT"

for gene in "${GENES[@]}"
do
  ASE_FILE="$gene/thoraxe/ases_table.csv"
  PATH_FILE="$gene/thoraxe/path_table.csv"

  if [[ ! -f "$ASE_FILE" || ! -f "$PATH_FILE" ]]; then
    echo "File not found: $ASE_FILE or $PATH_FILE, skipping $gene"
    continue
  fi

  counter1=0
  counter2=0
  counter3=0
  counter4=0

  # Get number of unique elements from column 5 of path_table.csv (skip header)
  total_paths=$(tail -n +2 "$PATH_FILE" | cut -d',' -f5 | sort | uniq | wc -l)

  if (( total_paths <= 0 )); then
    echo "No unique entries in column 5 of $PATH_FILE, skipping $gene"
    continue
  fi

  # Process ases_table.csv
  while IFS=, read -r col1 col2 col3 _ _ _ _ _ _ col10 _rest
  do
    if [[ "$col3" == "fully_alternative" ]]; then
      IFS='/' read -ra parts <<< "$col10"
      len=${#parts[@]}
      case $len in
        1) ((counter1++)) ;;
        2) ((counter2++)) ;;
        3) ((counter3++)) ;;
        4) ((counter4++)) ;;
        *) ;;  # ignore other cases
      esac
    fi
  done < <(tail -n +2 "$ASE_FILE")

  # Calculate percentages based on number of unique elements in column 5
  count1=$(echo "scale=2; $counter1*100/$total_paths" | bc)
  count2=$(echo "scale=2; $counter2*100/$total_paths" | bc)
  count3=$(echo "scale=2; $counter3*100/$total_paths" | bc)
  count4=$(echo "scale=2; $counter4*100/$total_paths" | bc)

  echo -e "${gene}\t${count1}\t${count2}\t${count3}\t${count4}" >> "$OUTPUT"
done
