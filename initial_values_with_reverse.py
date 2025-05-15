import os
import re
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

# Input files
gtf_file = "/home/abakarova/milica/IsoTools/Challenge2/iso_quant_PacBio_cDNA_IsoTools_final/human_simulation_PacBio_cDNA_long/models_ref_info.gtf"
genome_fasta = "/home/abakarova/milica/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
target_gene_id = "ENSG00000107643.16"

# Load the genome FASTA into a dictionary with chromosome names as keys
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# List to collect extracted exon data
records = []

# Counter to track exon numbers per transcript
exon_counter = defaultdict(int)

# Keywords to check in source field (2nd column)
keywords = ["HAVANA", "Bambu", "FLAIR", "PacBio_cDNA"]

# Parse GTF file
with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        if target_gene_id not in line:
            continue
        
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != "exon":
            continue
        
        chrom_raw = fields[0]
        chrom = chrom_raw.replace("chr", "")
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        source = fields[1]
        attributes = fields[8]
        
        attrs = dict(re.findall(r'(\S+) "([^"]+)"', attributes))
        
        gene_id = attrs.get("gene_id", "")
        transcript_id = attrs.get("transcript_id", "")
        exon_id = attrs.get("exon_id", "")

        # Find all tags present in the source field
        matched_tags = [tag for tag in keywords if tag in source]
        tag_suffix = "_".join(matched_tags) if matched_tags else ""

        if tag_suffix:
            gene_id += f"_{tag_suffix}"
            transcript_id += f"_{tag_suffix}"
            exon_id += f"_{tag_suffix}"

        # Increment and assign exon rank for this transcript
        exon_counter[transcript_id] += 1
        exon_rank = str(exon_counter[transcript_id])

        # Add transcript_id and exon rank to exon_id
        exon_id += f"_{transcript_id}_{exon_rank}"

        # Extract nucleotide sequence and adjust for strand
        if chrom in genome:
            seq = genome[chrom].seq[start - 1:end]
            if strand == "-":
                seq = seq.reverse_complement()
                start, end = end, start  # Swap start and end for consistency
            nucleotide_seq = str(seq)
        else:
            nucleotide_seq = "N/A"

        strand_value = 1 if strand == "+" else -1

        records.append([
            "homo_sapiens",
            gene_id,
            transcript_id,
            strand_value,
            exon_id,
            exon_rank,
            start,
            end,
            nucleotide_seq
        ])

# Create DataFrame
df = pd.DataFrame(records, columns=[
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "NucleotideSequence"
])

df["ExonRank"] = df["ExonRank"].astype(int)

# Separate positive and negative strands
df_pos = df[df["Strand"] == 1].sort_values(by=["TranscriptID", "ExonRank"])

df_neg = df[df["Strand"] == -1].sort_values(by=["TranscriptID", "ExonRank"], ascending=[True, False])

# Update exon ranks and IDs for negative strand

updated_neg = []
for transcript_id, group in df_neg.groupby("TranscriptID"):
    group = group.copy()
    group["ExonRank"] = range(1, len(group) + 1)
    group["ExonID"] = group.apply(
	lambda row: f"{row['GeneID']}_{row['TranscriptID']}_{row['ExonRank']}", axis=1
    )
    updated_neg.append(group)

if updated_neg:
    df_neg_updated = pd.concat(updated_neg)
else:
    # Kreiraj prazan DataFrame sa istim kolonama kao df_neg
    df_neg_updated = pd.DataFrame(columns=df_neg.columns)


# Combine all transcripts
df_final = pd.concat([df_pos, df_neg_updated])

# Final sort (optional)
df_final = df_final.sort_values(by=["TranscriptID", "ExonRank"])

# Save to CSV
output_csv = os.path.join(os.path.dirname(gtf_file), f"{target_gene_id}_exonsC.csv")
df_final.to_csv(output_csv, index=False)

print(f"Done. Table saved as: {output_csv}")
print(df_final.head())

