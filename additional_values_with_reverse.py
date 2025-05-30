import csv
from collections import defaultdict, namedtuple

ORF = namedtuple("ORF", ["start", "end", "frame"])

def find_orfs(seq, start_codons=["ATG"], stop_codons=["TAA", "TAG", "TGA"], min_length=30):
    seq = seq.upper()
    orfs = []
    for frame in range(3):
        start_pos = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if start_pos is None:
                if codon in start_codons:
                    start_pos = i
            else:
                if codon in stop_codons:
                    if (i + 3 - start_pos) >= min_length:
                        orfs.append(ORF(start=start_pos, end=i+3, frame=frame))
                    start_pos = None
    return orfs

def parse_csv(filename):
    transcripts = defaultdict(list)

    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            transcript_id = row['TranscriptID']
            exon = {
                'species': row['Species'],
                'gene_id': row['GeneID'],
                'transcript_id': transcript_id,
                'strand': row['Strand'],
                'exon_id': row['ExonID'],
                'exon_rank': int(row['ExonRank']),
                'start': int(row['ExonRegionStart']),
                'end': int(row['ExonRegionEnd']),
                'seq': row['NucleotideSequence'].strip().upper(),
                'cds_start': '.',
                'cds_end': '.',
                'start_phase': -1,
                'end_phase': -1,
            }
            transcripts[transcript_id].append(exon)

    for transcript_id, exons in transcripts.items():
        # descending sort
        exons.sort(key=lambda x: -x['exon_rank'])
        # rank update
        for new_rank, exon in enumerate(exons, 1):
            exon['exon_rank'] = new_rank
            # rotate start and end positions
            exon['start'], exon['end'] = exon['end'], exon['start']

        full_seq = ''.join([exon['seq'] for exon in exons])
        transcripts[transcript_id] = {
            'exons': exons,
            'sequence': full_seq,
        }

    return transcripts

def map_orf_to_cds(transcript, orf_start, orf_end):
    exons = transcript['exons']
    cumulative = 0

    for exon in exons:
        exon_len = len(exon['seq'])
        exon_start_in_transcript = cumulative
        exon_end_in_transcript = cumulative + exon_len

        overlap_start = max(orf_start, exon_start_in_transcript)
        overlap_end = min(orf_end, exon_end_in_transcript)

        if overlap_start < overlap_end:
            offset_start = overlap_start - exon_start_in_transcript
            offset_end = overlap_end - exon_start_in_transcript - 1 

            cds_start = exon['start'] + offset_start
            cds_end = exon['start'] + offset_end

            exon['cds_start'] = cds_start
            exon['cds_end'] = cds_end
        else:
            exon['cds_start'] = '.'
            exon['cds_end'] = '.'

        cumulative += exon_len

    return exons

def calculate_phase(transcript):
    exons = transcript['exons']
    prev_end_phase = None

    coding_exons = [exon for exon in exons if exon['cds_start'] != '.' and exon['cds_end'] != '.']
    if not coding_exons:
        return exons

    first_coding_exon = coding_exons[0]
    last_coding_exon = coding_exons[-1]

    for exon in exons:
        if exon['cds_start'] != '.' and exon['cds_end'] != '.':
            cds_start = exon['cds_start']
            cds_end = exon['cds_end']
            length = abs(cds_end - cds_start) + 1

            if exon is first_coding_exon:
                if cds_start == exon['start']:
                    exon['start_phase'] = 0
                else:
                    exon['start_phase'] = -1
                exon['end_phase'] = (exon['start_phase'] + length) % 3 if exon['start_phase'] != -1 else length % 3
                prev_end_phase = exon['end_phase']
            elif exon is not last_coding_exon:
                exon['start_phase'] = prev_end_phase
                exon['end_phase'] = (prev_end_phase + length) % 3
                prev_end_phase = exon['end_phase']
            if exon is last_coding_exon:
                if cds_end == exon['end']:
                    exon['end_phase'] = 0
                else:
                    exon['end_phase'] = -1
                if exon is not first_coding_exon:
                    exon['start_phase'] = prev_end_phase
        else:
            exon['start_phase'] = -1
            exon['end_phase'] = -1

    return exons

def write_csv(transcripts, output_file):
    with open(output_file, 'w', newline='') as f:
        fieldnames = [
            'Species', 'GeneID', 'TranscriptID', 'Strand',
            'ExonID', 'ExonRank', 'ExonRegionStart', 'ExonRegionEnd',
            'GenomicCodingStart', 'GenomicCodingEnd', 'StartPhase', 'EndPhase', 'NucleotideSequence'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for tid, transcript in transcripts.items():
            for exon in transcript['exons']:
                writer.writerow({
                    'Species': exon['species'],
                    'GeneID': exon['gene_id'],
                    'TranscriptID': exon['transcript_id'],
                    'Strand': exon['strand'],
                    'ExonID': exon['exon_id'],
                    'ExonRank': exon['exon_rank'],
                    'ExonRegionStart': exon['start'],
                    'ExonRegionEnd': exon['end'],
                    'GenomicCodingStart': exon['cds_start'],
                    'GenomicCodingEnd': exon['cds_end'],
                    'StartPhase': exon['start_phase'],
                    'EndPhase': exon['end_phase'],
                    'NucleotideSequence': exon['seq'],
                })

if __name__ == "__main__":
    input_csv = "/home/abakarova/milica/IsoTools/Challenge2/iso_quant_PacBio_cDNA_IsoTools_final/human_simulation_PacBio_cDNA_long/ENSG00000107643.16_exonsC.csv"
    output_csv = "/home/abakarova/milica/IsoTools/Challenge2/iso_quant_PacBio_cDNA_IsoTools_final/human_simulation_PacBio_cDNA_long/ENSG00000107643.16_exonsC_additional.csv"

    transcripts = parse_csv(input_csv)

    for tid, transcript in transcripts.items():
        orfs = find_orfs(transcript['sequence'])
        if orfs:
            longest_orf = max(orfs, key=lambda x: x.end - x.start)
            map_orf_to_cds(transcript, longest_orf.start, longest_orf.end)
            calculate_phase(transcript)

    write_csv(transcripts, output_csv)
