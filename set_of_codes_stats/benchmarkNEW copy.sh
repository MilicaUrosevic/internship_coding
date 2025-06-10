#!/bin/bash

GENES=(

"ENSMUSG00000019843"
#"ENSMUSG00000024985"

#"ENSMUSG00000018417"
#"ENSMUSG00000019647"
#"ENSMUSG00000020366"
#"ENSMUSG00000021936"
#"ENSMUSG00000022883"
#"ENSMUSG00000025020"
#"ENSMUSG00000027646"
#"ENSMUSG00000032366"
#"ENSMUSG00000052516"
#"ENSMUSG00000053825"
#"ENSMUSG00000057897"

)

SPECIES="mus_musculus,homo_sapiens"

# File for combined results
ALL_EVENTS="all_events_mouse.csv"
# Remove previous file if it exists
#rm -f "$ALL_EVENTS"

for GENE in "${GENES[@]}"; do
  echo "Processing gene: $GENE"

  transcript_query "$GENE" --specieslist $SPECIES --species "mus_musculus"
  ENS_DIR="${GENE}/Ensembl"

  if [[ -f "${ENS_DIR}/tsl.csv" ]]; then
    # Keep only homo_sapiens entries in tsl.csv
    awk -F',' 'NR==1 || $1 ~ /^mus_musculus$/' "${ENS_DIR}/tsl.csv" > "${ENS_DIR}/tsl_filtered.csv"
    mv "${ENS_DIR}/tsl_filtered.csv" "${ENS_DIR}/tsl.csv"
  else
    echo "Warning: ${ENS_DIR}/tsl.csv not found"
  fi

  if [[ -f "${ENS_DIR}/exonstable.tsv" ]]; then
    # Keep header and lines starting with ENSG (gene IDs)
    awk -F'\t' 'NR==1 || $1 ~ /^ENSMUS/' "${ENS_DIR}/exonstable.tsv" > "${ENS_DIR}/exonstable_filtered.tsv"
    mv "${ENS_DIR}/exonstable_filtered.tsv" "${ENS_DIR}/exonstable.tsv"
  else
    echo "Warning: ${ENS_DIR}/exonstable.tsv not found"
  fi

  if [[ -f "${ENS_DIR}/sequences.fasta" ]]; then
    # Filter out sequences whose header contains mus_musculus (remove mus_musculus entries)
    awk '
      /^>/ {
        if ($0 ~ /homo_sapiens/) {
          skip=1
        } else {
          skip=0
          print $0
        }
        next
      }
      skip==0 {print}
    ' "${ENS_DIR}/sequences.fasta" > "${ENS_DIR}/sequences_filtered.fasta"

    mv "${ENS_DIR}/sequences_filtered.fasta" "${ENS_DIR}/sequences.fasta"
  else
    echo "Warning: ${ENS_DIR}/sequences.fasta not found"
  fi

  python /Users/milicaurosevic/thoraxe/thoraxe/add_transcripts/add_transcripts.py "output/${GENE}_final.csv" "${ENS_DIR}"

  thoraxe -i "$GENE" -t 1 -s 5 --canonical_criteria "MinimumConservation,TranscriptLength"

  ASES_CSV="${GENE}/thoraxe/ases_table.csv"
  EVENTS_CSV="${GENE}/events_mouse.csv"
  PATH_TABLE="${GENE}/thoraxe/path_table.csv"

  if [[ ! -f "$ASES_CSV" ]]; then
    echo "Warning: $ASES_CSV not found, skipping events.csv creation"
    continue
  fi

  # === New code to determine longest_annotation flag ===
  longest_annotation_flag="F"
  if [[ -f "$PATH_TABLE" ]]; then
    # Get max value in col3 (numeric), ignoring header
    max_val=$(awk -F',' 'NR>1 {if($3>max) max=$3} END{print max}' "$PATH_TABLE")
    # Check if any row with col1 == GENE has col3 == max_val
    match=$(awk -F',' -v gene="$GENE" -v max="$max_val" 'NR>1 && $1==gene && $3==max {print "yes"}' "$PATH_TABLE")
    # Check if match is not empty (one or more matches)
    if [[ -n "$match" ]]; then
      longest_annotation_flag="T"
    fi
  else
    echo "Warning: $PATH_TABLE not found, setting longest_annotation=F"
  fi


  # Write header with new column longest_annotation
  echo -e "event\tflair\tbambu\trnabloom\tisotools\tannotation\tgene\tspecies\ttype\tlongest_annotation" > "$EVENTS_CSV"

  species="mus_musculus"

  rank=0
  tail -n +2 "$ASES_CSV" | while IFS=, read -r -a fields; do
    rank=$((rank + 1))

    canonical_field="${fields[8]}"
    alternative_field="${fields[9]}"

    function parse_tools {
      local field_str="$1"
      local flair="F"
      local bambu="F"
      local rnabloom="F"
      local isotools="F"
      local annotation="F"

      IFS='/' read -ra parts <<< "$field_str"
      for part in "${parts[@]}"; do
        part_upper=$(echo "$part" | tr '[:lower:]' '[:upper:]')

        if [[ "$part_upper" == *"FLAIR"* ]]; then flair="T"; fi
        if [[ "$part_upper" == *"BAMBU"* ]]; then bambu="T"; fi
        if [[ "$part_upper" == *"RNABLOOM"* ]]; then rnabloom="T"; fi
        if [[ "$part_upper" == *"PACBIO_CDNA"* ]]; then isotools="T"; fi

        if [[ "$part" =~ ^ENSMUSG[0-9]+$ ]]; then
          annotation="T"
        fi
      done

      echo -e "${flair}\t${bambu}\t${rnabloom}\t${isotools}\t${annotation}"
    }

    canonical_tools=$(parse_tools "$canonical_field")
    alternative_tools=$(parse_tools "$alternative_field")

    event_canonical="${GENE}_${rank}_canonical"
    event_alternative="${GENE}_${rank}_alternative"

    # Add longest_annotation_flag at the end of each row
    echo -e "${event_canonical}\t${canonical_tools}\t${GENE}\t${species}\tcanonical\t${longest_annotation_flag}"
    echo -e "${event_alternative}\t${alternative_tools}\t${GENE}\t${species}\talternative\t${longest_annotation_flag}"

  done >> "$EVENTS_CSV"

done

# Combine all events.csv into one file
echo "Combining all events.csv files into $ALL_EVENTS"

header_written=0
for GENE in "${GENES[@]}"; do
  EVENTS_CSV="${GENE}/events_mouse.csv"
  if [[ -f "$EVENTS_CSV" ]]; then
    if [[ $header_written -eq 0 ]]; then
      cat "$EVENTS_CSV" >> "$ALL_EVENTS"
      header_written=1
    else
      tail -n +2 "$EVENTS_CSV" >> "$ALL_EVENTS"
    fi
  fi
done

echo "All done. Combined events saved to $ALL_EVENTS"
