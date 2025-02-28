import re
import pandas as pd
import gffutils
import numpy as np

# Putanja do GTF fajla
gtf_file = "/Users/milicaurosevic/Desktop/TP53_test/all_zf_genes-vs-SRR8551562.gtf"
db_path = "gtf_database.db"

# Kreiranje baze podataka za GTF ako već ne postoji
try:
    db = gffutils.FeatureDB(db_path)
except:
    db = gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)

# Lista za skladištenje podataka
data = []

# Iteracija kroz egzone i vezivanje sa CDS regionima
for transcript in db.features_of_type("transcript"):
    transcript_id = transcript.attributes.get("transcript_id", [""])[0]
    gene_id = transcript.attributes.get("gene_id", [""])[0]
    strand = transcript.strand

    # Prikupimo sve CDS regione za dati transkript
    cds_regions = []
    for cds in db.children(transcript, featuretype="CDS", order_by="start"):
        cds_regions.append((cds.start, cds.end))

    # Sortiramo CDS regione po startnoj poziciji
    cds_regions.sort()

    # Iteriramo kroz egzone i dodeljujemo im CDS start/end vrednosti
    previous_end_phase = -1  # Početna vrednost za end_phase

    for exon in db.children(transcript, featuretype="exon", order_by="start"):
        exon_id = exon.id
        exon_rank = exon.attributes.get("exon_number", [""])[0]
        exon_start, exon_end = exon.start, exon.end

        # Pronalazimo CDS start i end unutar ovog egzona
        cds_start_coord, cds_end_coord = np.nan, np.nan
        for cds_start, cds_end in cds_regions:
            if exon_start <= cds_start <= cds_end <= exon_end:
                cds_start_coord, cds_end_coord = cds_start, cds_end
                break  # Uzimamo samo prvi koji upada u exon

        # Računanje start_phase i end_phase
        start_phase = -1
        end_phase = -1
        if not np.isnan(cds_start_coord) and not np.isnan(cds_end_coord):
            if previous_end_phase == -1:  # Prvi CDS
                start_phase = 0
                end_phase = (cds_end_coord - cds_start_coord + 1) % 3
            else:
                start_phase = previous_end_phase
                end_phase = (cds_end_coord - cds_start_coord + 1 + previous_end_phase) % 3

            previous_end_phase = end_phase  # Čuvamo za sledeći exon

        data.append([
            "Unknown", gene_id, transcript_id, strand, exon_id, exon_rank,
            exon_start, exon_end, cds_start_coord, cds_end_coord, start_phase, end_phase
        ])

# Pravljenje DataFrame-a
df = pd.DataFrame(data, columns=[
    "Species", "GeneID", "TranscriptID", "Strand", "ExonID", "ExonRank",
    "ExonRegionStart", "ExonRegionEnd", "GenomicCodingStart", "GenomicCodingEnd",
    "StartPhase", "EndPhase"
])

# Putanja za čuvanje CSV fajla
output_file = "/Users/milicaurosevic/Desktop/TP53_test/parsed_gtf_with_cds_fixed.csv"
df.to_csv(output_file, index=False)

print(f"Parsiranje završeno. Podaci su sačuvani u {output_file}")

# Prikaz prvih nekoliko redova da proveriš rezultat
print(df.head())
