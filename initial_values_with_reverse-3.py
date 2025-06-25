import os
import re
import sys
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

if len(sys.argv) != 2:
    print("Usage: python initial_values_with_reverse.py <target_gene_id>")
    sys.exit(1)

input_gene_id = sys.argv[1]
target_gene_id_short = input_gene_id.split(".")[0]  # Bez verzije


genome_fasta = "/home/abakarova/milica/human/benchmark_human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_files = [
    "/home/abakarova/milica/human/benchmark_human/Bambu.gtf",
    "/home/abakarova/milica/human/benchmark_human/IsoTools.gtf",
    "/home/abakarova/milica/human/benchmark_human/RNABloom.gtf",
    "/home/abakarova/milica/human/benchmark_human/FLAIR.gtf"
]

genome = {}
for record in SeqIO.parse(genome_fasta, "fasta"):
    chrom_id = record.id.split(".")[0]
    genome[chrom_id] = record

keywords = ["Bambu", "FLAIR", "PacBio_cDNA"]

records = []
exon_counter = defaultdict(int)

for gtf_file in gtf_files:
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            if line.startswith("#") or "\texon\t" not in line:
                continue

            match = re.search(r'gene_id "([^"]+)"', line)
            if not match:
                continue
            gtf_gene_id = match.group(1)
            gtf_gene_id_short = gtf_gene_id.split(".")[0]
            if gtf_gene_id_short != target_gene_id_short:
                continue

            fields = line.strip().split('\t')
            chrom_raw = fields[0]
            chrom = chrom_raw.replace("chr", "").split(".")[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            source = fields[1]
            attributes = fields[8]

            attrs = dict(re.findall(r'(\S+) "([^"]+)"', attributes))

            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            exon_id = attrs.get("exon_id", "")

            matched_tags = [tag for tag in keywords if tag in source]
            tag_suffix = "_".join(matched_tags) if matched_tags else "RNABLOOM"

            gene_id += f"_{tag_suffix}"
            transcript_id += f"_{tag_suffix}"
            exon_id += f"_{tag_suffix}"

            exon_counter[transcript_id] += 1
            exon_rank = str(exon_counter[transcript_id])
            exon_id += f"_{transcript_id}_{exon_rank}"

            if chrom in genome:
                seq = genome[chrom].seq[start - 1:end]
                if strand == "-":
                    seq = seq.reverse_complement()
                    start, end = end, start
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

if not records:
    print(f"No exons found for gene: {input_gene_id}")
    sys.exit(0)

df = pd.DataFrame(records, columns=[
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "NucleotideSequence"
])

df["ExonRank"] = df["ExonRank"].astype(int)

df_pos = df[df["Strand"] == 1].sort_values(by=["TranscriptID", "ExonRank"])
df_neg = df[df["Strand"] == -1].sort_values(by=["TranscriptID", "ExonRank"], ascending=[True, False])

updated_neg = []
for transcript_id, group in df_neg.groupby("TranscriptID"):
    group = group.copy()
    group["ExonRank"] = range(1, len(group) + 1)
    group["ExonID"] = group.apply(
        lambda row: f"{row['GeneID']}_{row['TranscriptID']}_{row['ExonRank']}", axis=1
    )
    updated_neg.append(group)

df_neg_updated = pd.concat(updated_neg) if updated_neg else pd.DataFrame(columns=df_neg.columns)

df_final = pd.concat([df_pos, df_neg_updated])
df_final = df_final.sort_values(by=["TranscriptID", "ExonRank"])

output_csv = os.path.join(os.getcwd(), f"{target_gene_id_short}_exons.csv")
df_final.to_csv(output_csv, index=False)

print(f"Done. Table saved as: {output_csv}")
print(df_final.head())

