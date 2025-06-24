# internship_coding

Repository for my M2 internship project:  
**Integrating long read data into Evolutionary Splicing Graphs (ESGs)**  
Supervisors: Elodie Laine (Sorbonne Université), Jean-Stéphane Varré (Université de Lille)

---

## Overview

This project explores how long-read transcript predictions can be integrated into Evolutionary Splicing Graphs using the ThorAxe framework.

Key goals:
- Evaluate long-read transcript prediction tools (FLAIR, Bambu, IsoTools, RNA-Bloom)
- Develop a module to integrate custom transcript annotations into ThorAxe
- Assess tool complementarity and their impact on splicing event discovery

---

## Main contributions

- `add_transcripts.py`: A Python module to add user-defined transcript structures to the ThorAxe input.
- Specification of a strict input CSV format derived from GTF and FASTA files.
- Comparative analysis of transcript-level and protein-level predictions from four tools.
- Integration of long-read–predicted isoforms into ESGs for human and mouse datasets.

---

## Usage

Run the module from the command line:

```bash
python add_transcripts.py path/to/transcripts.csv path/to/gene/Ensembl
