import gffutils
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

# Define input files
gtf_file = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.113-2.gtf"
genome_fasta = "/Users/milicaurosevic/Downloads/Taeniopygia_guttata.bTaeGut1_v1.p.dna.toplevel.fa"

# Define the GeneID to extract
target_gene_id = "ENSTGUG00000011517"

# Create a database from the GTF file (if not already created)
db_path = "ensembl.db"
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# Load the reference genome
genome_record = next(SeqIO.parse(genome_fasta, "fasta"))
genome_seq = genome_record.seq  # Extract the full sequence

# Prepare lists to store extracted data
data = []
cds_info = defaultdict(list)
utr5_lens = defaultdict(int)
utr3_lens = defaultdict(int)
cds_lengths = defaultdict(int)

# Preprocess all CDS, UTR5, and UTR3 regions
for feature in db.features_of_type(["CDS", "five_prime_utr", "three_prime_utr"]):
    transcript_id = feature.attributes.get("transcript_id", [""])[0]

    if feature.featuretype == "five_prime_utr":
        utr5_lens[transcript_id] += (feature.end - feature.start + 1)
    elif feature.featuretype == "three_prime_utr":
        utr3_lens[transcript_id] += (feature.end - feature.start + 1)
    elif feature.featuretype == "CDS":
        gene_id = feature.attributes.get("gene_id", [""])[0]
        if gene_id != target_gene_id:
            continue

        cds_start = feature.start
        cds_end = feature.end
        start_phase = int(feature.frame)  # Start phase from GTF

        # Compute end phase
        cds_length = cds_end - cds_start + 1
        end_phase = (start_phase + cds_length) % 3

        # Store CDS info sorted by start
        cds_info[transcript_id].append((cds_start, cds_end, start_phase, end_phase))
        cds_lengths[transcript_id] += cds_length

# Sort CDS regions per transcript to allow efficient searching
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

    # Retrieve matching CDS phase using sorted search
    matching_cds = None
    for cds_start, cds_end, start_phase, end_phase in cds_info.get(transcript_id, []):
        if cds_start >= exon_start and cds_end <= exon_end:
            matching_cds = (start_phase, end_phase)
            break  # No need to check further

    if matching_cds:
        start_phase, end_phase = matching_cds
    else:
        start_phase = "NA"
        end_phase = "NA"

    # Retrieve exon sequence efficiently
    exon_sequence = genome_seq[exon_start - 1 : exon_end]
    if strand == "-":
        exon_sequence = exon_sequence.reverse_complement()

    # Append to data list
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
        start_phase,
        end_phase,
        str(exon_sequence)
    ])

# Convert to DataFrame and save
columns = [
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd",
    "Start_Phase", "End_Phase", "NucleotideSequence"
]
df = pd.DataFrame(data, columns=columns)

output_file = f"optimized_exon_tablez_7_{target_gene_id}.csv"
df.to_csv(output_file, index=False)
print(f"Table saved as {output_file}")
