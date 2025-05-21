import pandas as pd

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

csv_file_path = '/Users/milicaurosevic/Downloads/MAPK8_exonsD_additional.csv'
df = pd.read_csv(csv_file_path)

df['TranscriptID'] = df['TranscriptID'].astype(str)

bambu = df['TranscriptID'].str.contains("Bambu", case=False)
flair = df['TranscriptID'].str.contains("FLAIR", case=False)
isotools = df['TranscriptID'].str.contains("PacBio_cDNA_models_ref_info", case=False)
rna_bloom = ~(bambu | flair | isotools)

groups = {
    "Bambu": df[bambu],
    "FLAIR": df[flairk],
    "IsoTools": df[isotools],
    "RNA Bloom": df[rna_bloom]
}

print("No of transcripts and proteins (different longest ORFs) per tool:\n")

for tool, group_df in groups.items():
    transcripts = group_df.groupby("TranscriptID")["NucleotideSequence"].apply(lambda x: ''.join(x.dropna())).to_dict()

    longest_orfs = set()
    for seq in transcripts.values():
        longest_orf = find_longest_orf(seq)
        if longest_orf:
            longest_orfs.add(longest_orf)

    print(f"{tool}:")
    print(f"  Number of transcripts: {len(transcripts)}")
    print(f"  Number of longest ORFs: {len(longest_orfs)}\n")
