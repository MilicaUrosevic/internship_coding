import gffutils
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Define input files
gtf_file = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.113-2.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"

# Define the GeneID to extract
target_gene_id = "ENSTGUG00000011517"  # Replace with the desired GeneID

# Create a database from the GTF file (if not already created)
db_path = "ensembl.db"
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# Load the reference genome (single chromosome)
genome_record = next(SeqIO.parse(genome_fasta, "fasta"))  # Fetch the only chromosome
genome_seq = genome_record.seq  # Extract the full sequence

# Prepare a list to store extracted exon data
data = []

# Store transcript UTR and CDS lengths
utr5_lens = defaultdict(int)
utr3_lens = defaultdict(int)
cds_lengths = defaultdict(int)

# Extract CDS, UTR5, and UTR3 regions
for feature in db.features_of_type(["CDS", "five_prime_utr", "three_prime_utr"]):
    transcript_id = feature.attributes.get("transcript_id", [""])[0]
    
    if feature.featuretype == "five_prime_utr":
        utr5_lens[transcript_id] += (feature.end - feature.start + 1)
    elif feature.featuretype == "three_prime_utr":
        utr3_lens[transcript_id] += (feature.end - feature.start + 1)
    elif feature.featuretype == "CDS":
        cds_lengths[transcript_id] += (feature.end - feature.start + 1)

# Extract exons only for the target GeneID
for exon in db.features_of_type("exon"):
    gene_id = exon.attributes.get("gene_id", [""])[0]
    
    if gene_id != target_gene_id:
        continue  # Skip if it’s not the target gene

    transcript_id = exon.attributes.get("transcript_id", [""])[0]
    exon_id = exon.id
    exon_rank = exon.attributes.get("exon_number", [""])[0]
    strand = exon.strand
    exon_start = exon.start
    exon_end = exon.end

    # Compute GenomicCodingStart and GenomicCodingEnd
    if transcript_id in cds_lengths and cds_lengths[transcript_id] > 0:
        if strand == "+":
            genomic_coding_start = exon_start + utr5_lens[transcript_id]
            genomic_coding_end = exon_end - utr3_lens[transcript_id]
        else:
            genomic_coding_start = exon_start + utr3_lens[transcript_id]
            genomic_coding_end = exon_end - utr5_lens[transcript_id]
    else:
        genomic_coding_start = "NA"
        genomic_coding_end = "NA"

    # Retrieve the exon sequence
    exon_sequence = genome_seq[exon_start - 1 : exon_end]  # Extract sequence

    if strand == "-":
        exon_sequence = exon_sequence.reverse_complement()

    # Append data to the list
    data.append([
        "taeniopygia_guttata",
        gene_id,
        transcript_id,
        strand,
        exon_id,
        exon_rank,
        exon_start,
        exon_end,
        genomic_coding_start,
        genomic_coding_end,
        str(exon_sequence)
    ])

# Convert to DataFrame
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd",
    "NucleotideSequence"
]
df = pd.DataFrame(data, columns=columns)

# Save results
output_file = f"exon_table_v5_{target_gene_id}.csv"
df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")
