#!/bin/bash
GENE_LIST="ENSG00000058404
ENSG00000197122
ENSG00000140416
ENSG00000139220
ENSG00000107643
ENSG00000187122
ENSG00000128641
ENSG00000148737
ENSG00000010810
ENSG00000169855
ENSG00000185008
ENSG00000092421
ENSG00000050748"

WORK_DIR="$(pwd)"
OUTPUT_DIR="$WORK_DIR/output"

mkdir -p "$OUTPUT_DIR"

echo "==> Pokrećem inicijalnu ekstrakciju egzona..."
for GENE_ID in $GENE_LIST; do
    echo "-> Obrada gena: $GENE_ID"
    python initial_values_with_reverse.py "$GENE_ID"
done

echo "==> Pokrećem dodatnu obradu i anotaciju ORF/CDS faza po fajlu..."

for CSV_FILE in "$WORK_DIR"/*_exons.csv; do
    BASENAME=$(basename "$CSV_FILE" _exons.csv)
    OUTPUT_CSV="$OUTPUT_DIR/${BASENAME}_final.csv"
    echo "-> Obrada fajla: $CSV_FILE"
    python additional_values_with_reverse.py "$CSV_FILE" "$OUTPUT_CSV"
done

echo "==> Sve završeno uspešno!"

