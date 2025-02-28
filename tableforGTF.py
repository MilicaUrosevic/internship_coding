import re
import pandas as pd

# path to gtf file
gtf_file = "/Users/milicaurosevic/Desktop/TP53_test/all_zf_genes-vs-SRR8551562.gtf"

# set the sp. name
species = "Unknown"

# list for saving data
data = []

# reading gtf file
with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue  # skip comments
        
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue  # skip if the line is not valid
        
        feature = fields[2]  # exon, CDS, etc.
        start = fields[3]
        end = fields[4]
        strand = fields[6]  # Strand (+ or -)
        
        # extract atteibutes from the last column
        attributes = fields[8]
        gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
        transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
        exon_id_match = re.search(r'exon_id "([^"]+)"', attributes)
        exon_rank_match = re.search(r'exon_number "([^"]+)"', attributes)
        
        if feature == "exon" and gene_id_match and transcript_id_match and exon_id_match and exon_rank_match:
            gene_id = gene_id_match.group(1)
            transcript_id = transcript_id_match.group(1)
            exon_id = exon_id_match.group(1)
            exon_rank = exon_rank_match.group(1)

            # Dodavanje podataka u listu
            data.append([species, gene_id, transcript_id, strand, exon_id, exon_rank, start, end])

# Pravljenje DataFrame-a
df = pd.DataFrame(data, columns=["Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank", "ExonRegionStart", "ExonRegionEnd"])

# Putanja za čuvanje CSV fajla
output_file = "/Users/milicaurosevic/Desktop/TP53_test/parsed_gtf_with_strand.csv"
df.to_csv(output_file, index=False)

print(f"Parsiranje završeno. Podaci su sačuvani u {output_file}")

# Prikaz prvih nekoliko redova ako želiš da proveriš
print(df.head())
