import gffutils
import pandas as pd
from Bio import SeqIO

#   Define input files
gtf_file = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.113-2.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"

#   Define the GeneID to extract
target_gene_id = "ENSTGUG00000027134"  # Replace with the desired GeneID

#   Create a database from the GTF file (if not already created)
db_path = "ensembl.db"
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

#  Load the reference genome, clean sequence IDs
genome_dict = {}
for record in SeqIO.parse(genome_fasta, "fasta"):
    clean_id = record.id.split()[0]  # Keep only the first part of the FASTA header
    genome_dict[clean_id] = record

# Prepare a list to store extracted exon data
data = []

# Extract exons only for the target GeneID
for exon in db.features_of_type("exon"):
    gene_id = exon.attributes.get("gene_id", [""])[0]
    
    # Filter only for the specified GeneID
    if gene_id != target_gene_id:
        continue  # Skip if it’s not the target gene

    transcript_id = exon.attributes.get("transcript_id", [""])[0]
    exon_id = exon.id
    exon_rank = exon.attributes.get("exon_number", [""])[0]
    strand = exon.strand
    exon_start = exon.start
    exon_end = exon.end

    # Retrieve the exon sequence
    seq_id = exon.seqid
    if seq_id in genome_dict:
        seq_record = genome_dict[seq_id]
        exon_sequence = seq_record.seq[exon_start - 1 : exon_end]  # Adjust for 0-based indexing
        if strand == "-":
            exon_sequence = exon_sequence.reverse_complement()
    else:
        exon_sequence = "NA"

    # Append data to the list
    data.append([
        "taeniopygia_guttata",  # Replace with the correct species name
        gene_id,
        transcript_id,
        strand,
        exon_id,
        exon_rank,
        exon_start,
        exon_end,
        str(exon_sequence)
    ])

# Convert to DataFrame
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "NucleotideSequence"
]
df = pd.DataFrame(data, columns=columns)

# Save results
output_file = f"exon_table_{target_gene_id}.csv"
df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")

