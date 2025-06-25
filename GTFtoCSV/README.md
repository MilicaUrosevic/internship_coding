# Transcript Table Generation Pipeline

This repository contains scripts to generate input CSV files used for the `add_transcripts` module. These CSV files are derived from GTF and FASTA files and include initial and additional transcript information.

## File Overview

- `makeTableS.sh` – main script that runs the pipeline
- `initial_values_with_reverse.py` – generates initial values from GTF and FASTA files
- `additional_values_with_reverse.py` – generates additional information for transcripts
- `GTFtoCSV/` – contains helper functions and conversion tools

## How to Run

The pipeline should be run in the following order:

1. Run the shell script `makeTableS.sh`  
   This script executes both the initial and additional value generation steps.

   ```bash
   bash makeTableS.sh
