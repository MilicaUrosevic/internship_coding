import re

# 1. Load the gene-to-chromosome mapping (without extensions)
gene_to_chr = {}

gene_to_chr_file = "/Users/milicaurosevic/Desktop/TP53_test/gictionary.txt"
bam_data_file = "/Users/milicaurosevic/Desktop/TP53_test/bam_data.txt"
output_file = "/Users/milicaurosevic/Desktop/TP53_test/absolute_splice_sites.txt"

# Load the dictionary while stripping spaces and removing extensions (.1, .2, etc.)
with open(gene_to_chr_file, "r") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            gene_id = parts[0].split('.')[0]  # Remove everything after the first "."
            chromosome = parts[1].strip()
            gene_to_chr[gene_id] = chromosome  # Store mapping

# 2. Function to extract absolute 3' splice site positions
def extract_absolute_splice_sites(start_pos, cigar_string):
    """Calculates absolute 3' splice site coordinates using the start position."""
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    genomic_pos = start_pos  # Initial genomic position of the read
    splice_sites = []  # List of 3' splice site positions
    
    for length, op in cigar_tuples:
        length = int(length)
        
        if op == 'M' or op == 'D':  
            # 'M' (match) and 'D' (deletion) move the genomic position forward
            genomic_pos += length
        
        elif op == 'N':  
            # 'N' represents an intron -> previous position is a 3' splice site
            splice_sites.append(genomic_pos)
            genomic_pos += length  # Skip over the intron

    return splice_sites

# 3. Read BAM data and map chromosomes using the gene-to-chromosome dictionary
with open(bam_data_file, 'r') as f, open(output_file, 'w') as out_f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 4:
            read_id = parts[0].strip()
            gene_id = parts[1].split('.')[0]  # Remove extension from gene ID
            start_pos = int(parts[2].strip())
            cigar = parts[3].strip()

            # Retrieve chromosome from the dictionary
            chromosome = gene_to_chr.get(gene_id, "Unknown")

            # Compute absolute 3' splice site positions
            splice_sites = extract_absolute_splice_sites(start_pos, cigar)

            # Save results
            out_f.write(f"{read_id} {chromosome} {splice_sites}\n")

print("Processing complete. Splice sites saved in", output_file)
