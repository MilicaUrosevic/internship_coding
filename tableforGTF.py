import re
import pandas as pd
import gffutils
import numpy as np
from Bio import SeqIO

# Path to the GTF file
gtf_file = "/Users/milicaurosevic/Desktop/TP53_test/all_zf_genes-vs-SRR8551562_with_chr.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"
db_path = "gtf_database.db"

# Load gene_id -> chromosome mapping
mapping_file = "/Users/milicaurosevic/Desktop/TP53_test/gictionary.txt"
map_dict = {}
with open(mapping_file, "r") as f:
    for line in f:
        fields = line.strip().split()
        if len(fields) >= 2:
            gene_id, chromosome = fields[0].split(".")[0], fields[1]  # Remove version from gene_id
            map_dict[gene_id] = chromosome

# Debug: Print first few entries of the dictionary
print("Sample of gene_id -> chromosome mapping:")
for i, (key, value) in enumerate(map_dict.items()):
    print(f"{key} -> {value}")
    if i >= 5:  # Print only first 5 for readability
        break

# Load genome sequence
genome_record = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Creating a database for the GTF file if it doesn't already exist
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# List to store data
data = []
# Iterating through transcripts and linking them with CDS regions
for transcript in db.features_of_type("transcript"):
    gene_id = transcript.chrom  # Move what was in the first column to gene_id
    transcript_id = transcript.attributes.get("transcript_id", [""])[0]
    original_gene_id = transcript.attributes.get("gene_id", [""])[0]  # Original incorrect gene_id
    chromosome = map_dict.get(gene_id.split(".")[0], "Unknown")  # Retrieve chromosome from the mapping

    # Debug: Print gene_id and its corresponding chromosome
    print(f"Processing gene_id: {gene_id}, Assigned chromosome: {chromosome}")
    
    strand = transcript.strand

    # Collect all CDS regions for the given transcript
    cds_regions = []
    for cds in db.children(transcript, featuretype="CDS", order_by="start"):
        cds_regions.append((cds.start, cds.end))
    
    # Sorting CDS regions by start position
    cds_regions.sort()

    # Iterating through exons and assigning them CDS start/end values
    previous_end_phase = -1  # Initial value for end_phase

    for exon in db.children(transcript, featuretype="exon", order_by="start"):
        exon_id = exon.id
        exon_rank = exon.attributes.get("exon_number", [""])[0]
        exon_start, exon_end = exon.start, exon.end  # Take start and end
        
        # Finding CDS start and end within this exon
        cds_start_coord, cds_end_coord = np.nan, np.nan
        for cds_start, cds_end in cds_regions:
            if exon_start <= cds_start <= cds_end <= exon_end:
                cds_start_coord, cds_end_coord = cds_start, cds_end
                break  # Take only the first one that fits within the exon

        # Calculating start_phase and end_phase
        start_phase = -1
        end_phase = -1
        if not np.isnan(cds_start_coord) and not np.isnan(cds_end_coord):
            if previous_end_phase == -1:  # First CDS
                start_phase = 0
                end_phase = (cds_end_coord - cds_start_coord + 1) % 3
            else:
                start_phase = previous_end_phase
                end_phase = (cds_end_coord - cds_start_coord + 1 + previous_end_phase) % 3

            previous_end_phase = end_phase  # Store for the next exon

        data.append([
            chromosome, gene_id, transcript_id, strand, exon_id, exon_rank,
            exon_start, exon_end, cds_start_coord, cds_end_coord, start_phase, end_phase 
        ])

# Creating a DataFrame
df = pd.DataFrame(data, columns=[
    "Chromosome", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd",
    "StartPhase", "EndPhase"
])

# Path to save the CSV file
output_file = "/Users/milicaurosevic/Desktop/TP53_test/parsed_gtf_corrected.csv"
df.to_csv(output_file, index=False)

print(f"Parsing completed. Data has been saved to {output_file}")

# Display the first few rows to check the result
print(df.head())
