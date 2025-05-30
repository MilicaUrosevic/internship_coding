import os
import re
import math
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import csv


gtf_files = {
    "bambu": "/home/abakarova/Downloads/Bambu/Challenge2/iso_quant_Bambu_cDNA_PacBio/mouse_simulation_cDNA_PacBio/models.gtf",
    "flair": "/home/abakarova/Downloads/Brooks/Challenge2/iso_quant_cdna_pacbio_brooks/mouse_sim_cdna_pacbio/models.gtf",
    "isotools": "/home/abakarova/milica/IsoTools/Challenge2/iso_quant_PacBio_cDNA_IsoTools_final/mouse_simulation_PacBio_cDNA_long/models_ref_info.gtf",
    "rnabloom": "/home/abakarova/Downloads/Birol/Challenge2/iso_quant_drna_ont_birol/mouse_simulation/lrgasp_gencode_vM27_sirvs.gtf",
}

genome_fasta = "/home/abakarova/milica/Mus_musculus.GRCm39.dna.toplevel.fa"
epsilon = 1


def parse_gtf_for_genes(gtf_path):
    gene_data = defaultdict(lambda: defaultdict(list))  # gene_id → transcript_id → list of exons
    strand_map = {}
    chrom_map = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            chrom, source, feature, start, end, _, strand, _, attributes = fields
            attrs = dict(re.findall(r'(\S+) "([^"]+)"', attributes))
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            start, end = int(start), int(end)

            gene_data[gene_id][transcript_id].append((start, end))
            strand_map[transcript_id] = strand
            chrom_map[transcript_id] = chrom.replace("chr", "")

    return gene_data, strand_map, chrom_map

# ----------------------
# longest orf function
# ----------------------

def find_longest_orf(seq, min_length=30):
    seq = seq.upper().replace("T", "U")
    start_codon = "AUG"
    stop_codons = {"UAA", "UAG", "UGA"}
    longest_orf = ""

    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(seq)-2, 3):
                    stop = seq[j:j+3]
                    if stop in stop_codons:
                        length = j + 3 - i
                        if length >= min_length:
                            orf_seq = seq[i:j+3]
                            if len(orf_seq) > len(longest_orf):
                                longest_orf = orf_seq
                        break
            i += 3
    return longest_orf

# ----------------------
# critical part starting here
# ----------------------

def are_orfs_equal(orf1, orf2):
    return orf1 == orf2

def count_unique_orfs(orf_list):
    unique_orfs = []
    for orf in orf_list:
        found = False
        for uorf in unique_orfs:
            if are_orfs_equal(orf, uorf):
                found = True
                break
        if not found:
            unique_orfs.append(orf)
    return len(unique_orfs)


def process_genes(genes_data):
    results = []
    for gene_id, (transcripts, strand_map, chrom_map, genome) in genes_data.items():
        transcript_count = len(transcripts)
        orf_sequences = []

        for transcript_id, exons in transcripts.items():
            exons = sorted(exons, key=lambda x: x[0])
            chrom = chrom_map.get(transcript_id, None)
            strand = strand_map.get(transcript_id, "+")
            if chrom is None or chrom not in genome:
                continue

            full_seq = ""
            for start, end in exons:
                seq = genome[chrom].seq[start-1:end]
                full_seq += str(seq)

            if strand == "-":
                full_seq = str(Seq(full_seq).reverse_complement())

            longest_orf = find_longest_orf(full_seq)  
            if longest_orf:
                orf_sequences.append(longest_orf)

        orf_count = count_unique_orfs(orf_sequences)  #here is the key problem
        results.append((gene_id, {"transcript_count": transcript_count, "orf_count": orf_count}))


    return results
    
# ----------------------
# critical part ends here
# ----------------------


# ----------------------
# parse gtf
# ----------------------

def parse_gtf_sequential(gtf_path, genome):
    print(f"Reading and parsing GTF file: {gtf_path}")
    gene_data, strand_map, chrom_map = parse_gtf_for_genes(gtf_path)

    print(f"Processing {len(gene_data)} genes sequentially.")
    results = []
    for gene_id, transcripts in gene_data.items():
        args = (gene_id, transcripts, strand_map, chrom_map, genome)
        res = process_genes({args[0]: args[1:]})[0]  
        results.append(res)

    result_dict = dict(results)
    return result_dict


# ----------------------
# main
# ----------------------

genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

tool_stats = {}
csv_rows = []

for tool, gtf_path in gtf_files.items():
    print(f"Parsing {tool} GTF sequentially...")
    tool_stats[tool] = parse_gtf_sequential(gtf_path, genome)

    for gene_id, stats in tool_stats[tool].items():
        transcript_count = stats["transcript_count"]
        orf_count = stats["orf_count"]

        for _ in range(transcript_count):
            csv_rows.append({
                "transcript_or_protein_count": 1,
                "gene_id": gene_id,
                "method": tool,
                "species": "mus_musculus",
                "polymer": "rna"
            })

        for _ in range(orf_count):
            csv_rows.append({
                "transcript_or_protein_count": 1,
                "gene_id": gene_id,
                "method": tool,
                "species": "mus_musculus",
                "polymer": "protein"
            })

# ----------------------
# compare with bambu
# ----------------------

data_combined = []

for gene_id in tool_stats["bambu"]:
    bambu_transcripts = tool_stats["bambu"].get(gene_id, {}).get("transcript_count", 0)
    bambu_orfs = tool_stats["bambu"].get(gene_id, {}).get("orf_count", 0)

    for tool in gtf_files:
        if tool == "bambu":
            continue
        tool_transcripts = tool_stats[tool].get(gene_id, {}).get("transcript_count", 0)
        tool_orfs = tool_stats[tool].get(gene_id, {}).get("orf_count", 0)

        log2_trans = math.log2((bambu_transcripts + epsilon) / (tool_transcripts + epsilon))
        log2_orfs = math.log2((bambu_orfs + epsilon) / (tool_orfs + epsilon))

        data_combined.append({"gene_id": gene_id, "tool": tool, "log2_ratio": log2_trans, "type": "transcript"})
        data_combined.append({"gene_id": gene_id, "tool": tool, "log2_ratio": log2_orfs, "type": "orf"})

# ----------------------
# plot boxplot
# ----------------------

def plot_combined(data, filename):
    df = pd.DataFrame(data)

    # Mapa sa "prikaznim" imenima za alat (koriste se na x osi)
    display_names = {
        "bambu": "Bambu",
        "flair": "FLAIR",
        "isotools": "IsoTools",
        "rnabloom": "RNA Bloom"
    }

    # Dodaj novu kolonu sa prikaznim imenima alata za lakšu upotrebu u plotu
    df['tool_display'] = df['tool'].map(display_names)

    plt.figure(figsize=(10, 7))
    
    # Definiši paletu po tipu ("type")
    custom_palette = {"transcript": "lightblue",  
                      "orf": "mistyrose"}        

    sns.boxplot(x="tool_display", y="log2_ratio", hue="type", data=df, palette=custom_palette)
    
    plt.axhline(0, linestyle="--", color="black")
    plt.title("log2(Bambu / Tool): Transcript and ORF counts per gene for Mus musculus")
    plt.ylabel("log2 Ratio")
    plt.xlabel("Tool (compared to Bambu)")
    plt.legend(title="Data type")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved combined plot: {filename}")

# Poziv funkcije:
plot_combined(data_combined, "log2_transcript_orf_counts_per_gene_human_aaaaa.png")




def plot_log2_total_polymers(tool_stats, epsilon=1, filename="log2_total_polymers_per_tool_aaaaa.png"):
    # Mapa sa "prikaznim" imenima za alat (koriste se na x osi)
    display_names = {
        "bambu": "Bambu",
        "flair": "FLAIR",
        "isotools": "IsoTools",
        "rnabloom": "RNA Bloom"
    }

    # Lista sa internim imenima alata (mala slova), redosled je bitan za x osu i sumiranje
    tools = ["bambu", "flair", "isotools", "rnabloom"]

    total_transcripts = []
    total_orfs = []

    for tool in tools:
        # Sumiranje vrednosti za svaki alat koristeći mala imena
        sum_transcripts = sum(stats["transcript_count"] for stats in tool_stats[tool].values())
        sum_orfs = sum(stats["orf_count"] for stats in tool_stats[tool].values())

        # log2 transformacija
        total_transcripts.append(math.log2(sum_transcripts + epsilon))
        total_orfs.append(math.log2(sum_orfs + epsilon))

    x = range(len(tools))

    plt.figure(figsize=(10, 6))
    plt.bar(x, total_transcripts, width=0.4, label="log2 Total Transcripts", align='center', color='lightblue')
    plt.bar([i + 0.4 for i in x], total_orfs, width=0.4, label="log2 Total ORFs", align='center', color='mistyrose')

    # Koristimo prikazna imena za x oznake
    plt.xticks([i + 0.2 for i in x], [display_names[t] for t in tools], rotation=45, ha='right')

    plt.ylabel("log2 total polymer counts")
    plt.xlabel("Tools")
    plt.title(r"log2 total polymer counts obtained by different tools for Mus musculus", fontsize=14)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

    print(f"Saved bar plot of total log2 counts per tool as {filename}")

# Poziv funkcije:
plot_log2_total_polymers(tool_stats)


# ----------------------
# csv
# ----------------------

with open("transcript_protein_counts_human_orf.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["transcript_or_protein_count", "gene_id", "method", "species", "polymer"])
    writer.writeheader()
    writer.writerows(csv_rows)

