import gffutils
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Define input files
gtf_file = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.113-2.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"
target_gene_id = "ENSTGUG00000008103"

db_path = "ensembl.db"
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# Load the reference genome
genome_record = next(SeqIO.parse(genome_fasta, "fasta"))
genome_seq = genome_record.seq  # Extract the full sequence

cds_info = defaultdict(list)
exon_data = []

# Collect all CDS information for the target gene
for feature in db.features_of_type("CDS"):
    if feature.attributes.get("gene_id", [""])[0] == target_gene_id:
        transcript_id = feature.attributes.get("transcript_id", [""])[0]
        cds_info[transcript_id].append((feature.start, feature.end))

# Sort CDS regions per transcript
for transcript in cds_info:
    cds_info[transcript].sort()

# Extract exons only for the target gene
for exon in db.features_of_type("exon"):
    if exon.attributes.get("gene_id", [""])[0] != target_gene_id:
        continue
    
    transcript_id = exon.attributes.get("transcript_id", [""])[0]
    exon_rank = exon.attributes.get("exon_number", [""])[0]
    strand = exon.strand
    exon_start, exon_end = exon.start, exon.end
    
    # Find the first CDS within this exon
    cds_start_coord, cds_end_coord = "NA", "NA"
    if transcript_id in cds_info:
        for cds_start, cds_end in cds_info[transcript_id]:
            if exon_start <= cds_start <= cds_end <= exon_end:
                cds_start_coord, cds_end_coord = cds_start, cds_end
                break
    
    # Retrieve exon sequence
    exon_sequence = genome_seq[exon_start - 1 : exon_end]
    if strand == "-":
        exon_sequence = exon_sequence.reverse_complement()
    
    exon_data.append([
        "taeniopygia_guttata", target_gene_id, transcript_id, strand, exon.id, exon_rank,
        exon_start, exon_end, cds_start_coord, cds_end_coord, str(exon_sequence)
    ])

# Compute start and end phase
processed_data = []
previous_end_phase = -1
start_phase = -1
for row in exon_data:
    species, gene_id, transcript_id, strand, exon_id, exon_rank, exon_start, exon_end, cds_start, cds_end, seq = row
    end_phase = -1
    
    if cds_start != "NA" and cds_end != "NA":
        if previous_end_phase == -1:
            start_phase = -1
            end_phase = (int(cds_end) - int(cds_start) + 1) % 3
            previous_end_phase = end_phase
        else:
            start_phase = previous_end_phase
            end_phase = (int(cds_end) - int(cds_start) + 1 + previous_end_phase) % 3
            previous_end_phase = end_phase
    processed_data.append(row + [start_phase, end_phase])

# Convert to DataFrame and save
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd", 
    "NucleotideSequence", "StartPhase", "EndPhase"
]
df = pd.DataFrame(processed_data, columns=columns)

output_file = f"optimized_exon_table_v8_{target_gene_id}.csv"
df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")
