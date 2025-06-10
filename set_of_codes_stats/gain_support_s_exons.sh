#!/bin/bash

# Define the list of genes
GENES=(
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

# Define suffixes of interest
SUFFIXES="_Bambu|_FLAIR|_cDNA|_RNABLOOM"

# Temporary directory to store intermediate results
TMP_DIR=$(mktemp -d)

# Process each gene's s_exon_table.csv file
for gene in "${GENES[@]}"; do
  FILE="${gene}/thoraxe/s_exon_table.csv"
  
  # If file doesn't exist, write missing values
  if [[ ! -f "$FILE" ]]; then
    echo "missing,missing,missing,missing" > "$TMP_DIR/${gene}_support"
    echo "missing,missing,missing,missing" > "$TMP_DIR/${gene}_gain"
    continue
  fi

  # Process the file with AWK
  gawk -v gene="$gene" -v SUFFIXES="$SUFFIXES" '
  BEGIN { FS=OFS="," }
  {
    s_exon=$6;
    geneID=$2;
    data[s_exon][geneID]=1;
    genes[s_exon]=1;
  }
  END {
    for (i = 1; i <= 4; i++) {
      support[i] = 0;
      gain[i] = 0;
    }
    total = 0;

    for (s in genes) {
      split("", gids); suf_count=0; has_gene=0;

      # Check which GeneIDs are associated with the current s_exon
      for (g in data[s]) {
        gids[g] = 1;
        if (g == gene) {
          has_gene = 1;
        } else if (g ~ SUFFIXES) {
          suf_count++;
        }
      }

      total++;

      # Case 1: gene not in GeneIDs, but suffix count = 1, 2, 3, or 4 (support)
      if (!has_gene && suf_count >= 1 && suf_count <= 4) {
        support[suf_count]++;
      }

      # Case 2: gene is in GeneIDs, and suffix count = 1, 2, 3, or 4 (gain)
      if (has_gene && suf_count >= 1 && suf_count <= 4) {
        gain[suf_count]++;
      }
    }

    # Normalize results by number of distinct s_exons and scale by 6
    for (i = 1; i <= 4; i++) {
      if (total > 0) {
        printf "%.4f%s", support[i]/total*100, (i<4?",":"\n")
      } else {
        printf "0.0000%s", (i<4?",":"\n")
      }
    }

    for (i = 1; i <= 4; i++) {
      if (total > 0) {
        printf "%.4f%s", gain[i]/total*100, (i<4?",":"\n")
      } else {
        printf "0.0000%s", (i<4?",":"\n")
      }
    }
  }
  ' "$FILE" | {
    read support_line
    read gain_line
    echo "$support_line" > "$TMP_DIR/${gene}_support"
    echo "$gain_line" > "$TMP_DIR/${gene}_gain"
  }
done

# Create CSV headers
echo "Gene,1,2,3,4" > support_matrix_inverted.csv
echo "Gene,1,2,3,4" > gain_matrix_inverted.csv

# Write results from temp files into final CSVs
for gene in "${GENES[@]}"; do
  support_values=$(cat "$TMP_DIR/${gene}_support" 2>/dev/null || echo "missing,missing,missing,missing")
  gain_values=$(cat "$TMP_DIR/${gene}_gain" 2>/dev/null || echo "missing,missing,missing,missing")
  echo "$gene,$support_values" >> support_matrix_inverted.csv
  echo "$gene,$gain_values" >> gain_matrix_inverted.csv
done

# Clean up temporary files
rm -r "$TMP_DIR"
