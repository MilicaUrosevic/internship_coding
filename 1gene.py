import gffutils
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Define input files
gtf_file = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.113-2.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"

target_gene_id = "ENSTGUG00000011517"

db_path = "ensembl.db"
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# Load the reference genome
genome_record = next(SeqIO.parse(genome_fasta, "fasta"))
genome_seq = genome_record.seq  # Extract the full sequence

data = []
cds_info = defaultdict(list)

# Collect all CDS information for the target gene
for feature in db.features_of_type("CDS"):
    gene_id = feature.attributes.get("gene_id", [""])[0]
    if gene_id != target_gene_id:
        continue
    
    transcript_id = feature.attributes.get("transcript_id", [""])[0]
    cds_info[transcript_id].append((feature.start, feature.end))

# Sort CDS regions per transcript
for transcript in cds_info:
    cds_info[transcript].sort()

# Extract exons only for the target GeneID
for exon in db.features_of_type("exon"):
    gene_id = exon.attributes.get("gene_id", [""])[0]
    if gene_id != target_gene_id:
        continue
    
    transcript_id = exon.attributes.get("transcript_id", [""])[0]
    exon_id = exon.id
    exon_rank = exon.attributes.get("exon_number", [""])[0]
    strand = exon.strand
    exon_start = exon.start
    exon_end = exon.end
    
    # Find the first CDS within this exon and before the next one
    cds_start_coord = "NA"
    cds_end_coord = "NA"
    
    if transcript_id in cds_info:
        for cds_start, cds_end in cds_info[transcript_id]:
            if cds_start >= exon_start and cds_end <= exon_end:
                cds_start_coord = cds_start
                cds_end_coord = cds_end
                break
    
    # Retrieve exon sequence
    exon_sequence = genome_seq[exon_start - 1 : exon_end]
    if strand == "-":
        exon_sequence = exon_sequence.reverse_complement()
    
    # Append to data list
    data.append([
        "taeniopygia_guttata", gene_id, transcript_id, strand, exon_id, exon_rank,
        exon_start, exon_end, cds_start_coord, cds_end_coord, str(exon_sequence)
    ])

# Convert to DataFrame and save
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd", "NucleotideSequence"
]
df = pd.DataFrame(data, columns=columns)

output_file = f"optimized_exon_table_v8_{target_gene_id}.csv"
df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")
