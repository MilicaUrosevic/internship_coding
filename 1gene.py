import gffutils
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import numpy as np

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
genome_record = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))


cds_info = defaultdict(list)
stop_codon_info = defaultdict(list)
exon_data = []

for feature in db.features_of_type("CDS"):
    if feature.attributes.get("gene_id", [""])[0] == target_gene_id:
        transcript_id = feature.attributes.get("transcript_id", [""])[0]
        cds_info[transcript_id].append((feature.chrom, feature.start, feature.end))


# Extract stop codons for the target gene
for feature in db.features_of_type("stop_codon"):
    if feature.attributes.get("gene_id", [""])[0] == target_gene_id:
        transcript_id = feature.attributes.get("transcript_id", [""])[0]
        stop_codon_info[transcript_id].append((feature.start, feature.end))


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
    exon_start, exon_end, chrom = exon.start, exon.end, exon.chrom
    
    # Find the first CDS within this exon
    cds_start_coord, cds_end_coord = np.nan, np.nan
    if transcript_id in cds_info:
        for chrom, cds_start, cds_end in cds_info[transcript_id]:
            if exon_start <= cds_start <= cds_end <= exon_end:
                cds_start_coord, cds_end_coord = cds_start, cds_end
                # Check if a stop codon follows this CDS
                if transcript_id in stop_codon_info:
                    for stop_codon_start, stop_codon_end in stop_codon_info[transcript_id]:
                        if cds_end + 1 == stop_codon_start:
                            cds_end_coord = stop_codon_end  # Extend CDS end to stop codon end
                break

    
    # Retrieve exon sequence
    genome_seq = genome_record[chrom].seq
    exon_sequence = genome_seq[exon_start - 1 : exon_end]
    if strand == "-":
        exon_sequence = exon_sequence.reverse_complement()

    if strand == "-":
        strand = -1
    else:
        strand = 1
    
    exon_data.append([
        "taeniopygia_guttata", target_gene_id, transcript_id, strand, exon.id, exon_rank,
        exon_start, exon_end, cds_start_coord, cds_end_coord, str(exon_sequence)
    ])


# Calculate start and end phase
processed_data = []
previous_end_phase = -1  # initialization for previous end phase

for i, row in enumerate(exon_data):
    species, gene_id, transcript_id, strand, exon_id, exon_rank, exon_start, exon_end, cds_start, cds_end, seq = row
    start_phase = end_phase = -1  

    #check if it is the last exon of a transcript
    is_last_exon = (i == len(exon_data) - 1) or (exon_data[i + 1][2] != transcript_id)
    
    if cds_start is not None and cds_end is not None and not np.isnan(cds_start) and not np.isnan(cds_end):
        if previous_end_phase == -1:
            start_phase = previous_end_phase
            end_phase = (int(cds_end) - int(cds_start) + 1) % 3
        else:
            start_phase = previous_end_phase
            end_phase = (int(cds_end) - int(cds_start) + 1 + previous_end_phase) % 3

        if is_last_exon:
            end_phase = -1
            if cds_end == exon_end:
                end_phase = 0

        previous_end_phase = end_phase
    else:
        previous_end_phase = -1

    new_row = list(row)  
    new_row.insert(-1, start_phase) 
    new_row.insert(-1, end_phase)    
    
    processed_data.append(new_row)  

# Convert to DataFrame and save
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd", 
    "StartPhase", "EndPhase", "NucleotideSequence"
]
df = pd.DataFrame(processed_data, columns=columns)

output_file = f"zebrafinch6{target_gene_id}.csv"
df.to_csv(output_file, index=False, sep=",")
print(f"Table saved as {output_file}")