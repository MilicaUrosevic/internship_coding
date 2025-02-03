"""
Functions to add user-defined transcript data from a CSV file to the
previously downloaded Ensembl data.
"""

import argparse
import collections
import datetime
import logging
import os

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet

from tabulate import tabulate
from thoraxe.transcript_info import read_exon_file
from thoraxe.version import __version__


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse add_transcripts' command line arguments and
    returns the argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="add_transcripts",
        description="""
        add_transcripts downloads add user-defined transcripts to
        the data previously download from Ensembl.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input',
                        help='Input CSV containing the transcript data.',
                        type=str)
    parser.add_argument(
        'ensembl',
        help='Path to the previously download Ensembl data for a gene, '
        'i.e. the path to the Ensembl directory created by '
        'transcript_query',
        type=str)
    parser.add_argument('-v',
                        '--verbose',
                        help='Print detailed progress.',
                        action='store_true')
    parser.add_argument('--version',
                        action='version',
                        version=f'ThorAxe version {__version__}')

    return parser


COLUMN_NAMES = {
    'Species',
    'GeneID',
    'TranscriptID',
    'Strand',
    'ExonID',
    'ExonRank',
    'ExonRegionStart',
    'ExonRegionEnd',
    'GenomicCodingStart',
    'GenomicCodingEnd',
    'StartPhase',
    'EndPhase',
    'NucleotideSequence',
}

Paths = collections.namedtuple('Paths',
                               ['input', 'tsl', 'exontable', 'sequences'])


def get_output_paths(args):
    """
    Return a `Paths` `namedtuple` containing the input and output paths for
    `add_transcripts`.
    """
    return Paths(args.input, os.path.join(args.ensembl, 'tsl.csv'),
                 os.path.join(args.ensembl, 'exonstable.tsv'),
                 os.path.join(args.ensembl, 'sequences.fasta'))


def read_fasta(sequences_path):
    """
    Returns a list of `SeqRecord`s for the fasta file at `sequences_path`.
    """
    return list(SeqIO.parse(sequences_path, "fasta"))


def _write_fasta(sequences_path, seqrecords):
    SeqIO.write(seqrecords, sequences_path, "fasta")


def _get_exon_sequences(seqrecords):
    return {
        seqrecord.description.split(' ')[1]: str(seqrecord.seq)
        for seqrecord in seqrecords
    }


def _in(value, column):
    return any(value == column)


def _check_col_names(input_df):
    for col in COLUMN_NAMES:
        if col not in input_df.columns:
            raise Exception(f'The {col} column is missing in the input table.')
    for col in input_df.columns:
        if col not in COLUMN_NAMES:
            logging.warning('The %s column is not going to be used', col)


def _check_transcript(input_df, exontable):
    for transcript in input_df['TranscriptID'].unique():
        if _in(transcript, exontable['TranscriptID']):
            raise Exception(
                f'TranscriptID {transcript} is already in the Ensembl data!')


def _format_df(data):
    return tabulate(data,
                    headers='keys',
                    tablefmt='psql',
                    floatfmt=".0f",
                    showindex=False)


def _check_exon(input_df, exontable):
    exon_columns = [
        'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart',
        'GenomicCodingEnd', 'StartPhase', 'EndPhase'
    ]
    for exon in input_df['ExonID'].unique():
        if _in(exon, exontable['ExonID']):
            logging.warning(
                "%s is on the Ensembl data; the new sequence won't be used.",
                exon)
            previous_exon = exontable[exontable.ExonID ==
                                      exon][exon_columns].drop_duplicates()
            new_exon = input_df[input_df.ExonID ==
                                exon][exon_columns].drop_duplicates()
            if not (previous_exon.values == new_exon.values).all():
                raise Exception(f'''
Exon {exon} is already in the Ensembl data with different values!
                    
Exon data at Ensembl: 
{_format_df(previous_exon)}

New exon data:
{_format_df(new_exon)}''')


def _check_species_name(input_df):
    for row in input_df.itertuples():
        example = 'e.g. homo_sapiens.'
        name = row.Species
        if not name.islower():
            raise Exception(
                f'Error with {name}: Species name should be lowercase, {example}'
            )
        if len(name.split('_')) != 2:
            raise Exception(
                f'Error with {name}:' +
                'Species name should be binomial, and the terms should be '
                f'separated by underscore, {example}')


def _check_pair_order(input_df):
    for row in input_df.itertuples():
        if row.ExonRegionStart > row.ExonRegionEnd:
            raise Exception(
                'ExonRegionStart should be lower than ExonRegionEnd.')
        if row.GenomicCodingStart > row.GenomicCodingEnd:
            raise Exception(
                'GenomicCodingStart should be lower than GenomicCodingEnd.')
        if len(row.NucleotideSequence) != (row.ExonRegionEnd -
                                           row.ExonRegionStart + 1):
            raise Exception(
                "The sequence length estimated using ExonRegionStart and "
                "ExonRegionEnd disagrees with the actual NucleotideSequence "
                "length.")


def check_input(input_df, exontable):
    """
    It raises an error if there is a problem with the user-defined transcripts.
    """
    _check_col_names(input_df)
    _check_species_name(input_df)
    _check_transcript(input_df, exontable)
    _check_pair_order(input_df)
    _check_exon(input_df, exontable)

"""
def add_to_exontable(input_df, exontable):
    """"""
    It adds the new exons and transcripts to the Ensembl's exontable,
    computing cDNA start and end values based on the given conditions.
    """"""
    
    # Sort input dataframe by TranscriptID and ExonRank
    input_df = input_df.sort_values(by=['TranscriptID', 'ExonRank']).reset_index(drop=True)
    
    transcript_id = None
    exon_len_accumulator = 0  # Accumulates exon lengths when genomic_len is 0
    
    for row in input_df.itertuples():
        genomic_len = row.GenomicCodingEnd - row.GenomicCodingStart
        exon_len = row.ExonRegionEnd - row.ExonRegionStart
        
        if transcript_id != row.TranscriptID:
            transcript_id = row.TranscriptID
            exon_len_accumulator = 0  # Reset for new transcript
        
        if genomic_len == 0:
            cdna_start = 0
            cdna_end = 0
            exon_len_accumulator += exon_len  # Accumulate exon length
        else:
            cdna_start = exon_len_accumulator + abs(row.GenomicCodingStart - row.ExonRegionStart) + 1
            cdna_end = cdna_start + genomic_len
            exon_len_accumulator += exon_len  # Update exon length accumulator
        
        exontable = exontable.append(
            {
                'GeneID': row.GeneID,
                'TranscriptID': row.TranscriptID,
                'ProteinID': f"{row.TranscriptID}_PROTEIN",
                'Strand': row.Strand,
                'ExonID': row.ExonID,
                'ExonRegionStart': row.ExonRegionStart,
                'ExonRegionEnd': row.ExonRegionEnd,
                'ExonRank': row.ExonRank,
                'cDNA_CodingStart': cdna_start,
                'cDNA_CodingEnd': cdna_end,
                'GenomicCodingStart': row.GenomicCodingStart,
                'GenomicCodingEnd': row.GenomicCodingEnd,
                'StartPhase': row.StartPhase,
                'EndPhase': row.EndPhase,
            },
            ignore_index=True
        )
    
    return exontable"""

def calculate_cDNA_start_end(df):
    # Sort by TranscriptID and ExonRank - done
    df = df.sort_values(by=['TranscriptID', 'ExonRank']).copy()
    
    # Convert necessary columns to numeric type
    numeric_columns = ['ExonRank', 'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart', 'GenomicCodingEnd', 'Strand']
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Initialize new columns
    df['cDNA_CodingStart'] = np.nan
    df['cDNA_CodingEnd'] = np.nan
    
    # Iterate through groups by TranscriptID
    for transcript_id, transcript_group in df.groupby('TranscriptID'):
        transcript_group = transcript_group.sort_values(by='ExonRank')
        
        cDNA_position = 1  # Reset start position for each transcript
        cumulative_exon_length = 0  # Accumulate exon length before the first coding exon
        found_coding_exon = False  # Track if the first coding exon is encountered
        
        for idx, row in transcript_group.iterrows():
            if not found_coding_exon:
                if pd.notna(row['GenomicCodingStart']) and pd.notna(row['GenomicCodingEnd']):
                    # Adjust calculation based on strand direction 
                    coding_length = abs(row['GenomicCodingEnd'] - row['GenomicCodingStart'])
                   
                    if row['Strand'] == 1:
                        df.at[idx, 'cDNA_CodingStart'] = cumulative_exon_length + abs(row['GenomicCodingStart'] - row['ExonRegionStart']) + row['ExonRank']
                        df.at[idx, 'cDNA_CodingEnd'] = df.at[idx, 'cDNA_CodingStart'] + coding_length
                        cDNA_position = df.at[idx, 'cDNA_CodingEnd'] + 1
                        found_coding_exon = True  # The first coding exon has been found
                    elif row['Strand'] == -1:
                        df.at[idx, 'cDNA_CodingStart'] = cumulative_exon_length + abs(row['GenomicCodingEnd'] - row['ExonRegionEnd']) + row['ExonRank']
                        df.at[idx, 'cDNA_CodingEnd'] = df.at[idx, 'cDNA_CodingStart'] + coding_length
                        cDNA_position = df.at[idx, 'cDNA_CodingEnd'] + 1
                        found_coding_exon = True  # The first coding exon has been found
                else:
                    # Accumulate exon length before the coding exon
                    cumulative_exon_length += abs(row['ExonRegionEnd'] - row['ExonRegionStart'])
            else:
                # Adjust calculation based on strand direction
                
                coding_length = abs(row['GenomicCodingEnd'] - row['GenomicCodingStart'])
                df.at[idx, 'cDNA_CodingStart'] = cDNA_position
                df.at[idx, 'cDNA_CodingEnd'] = cDNA_position + coding_length
                cDNA_position = df.at[idx, 'cDNA_CodingEnd'] + 1
    
    return df


def add_to_exontable(input_df, exontable):
    # First, calculate cDNA start and end values
    input_df = calculate_cDNA_start_end(input_df)
    
    # Remove NaN values and convert numeric columns
    input_df.dropna(inplace=True)
    numeric_columns = [
        'ExonRegionStart', 'ExonRegionEnd', 'GenomicCodingStart', 
        'GenomicCodingEnd', 'StartPhase', 'EndPhase', 'cDNA_CodingStart', 'cDNA_CodingEnd'
    ]
    
    for col in numeric_columns:
        if col in input_df.columns:
            input_df[col] = pd.to_numeric(input_df[col], errors='coerce').astype('Int64')
    
    # Sort data
    input_df = input_df.sort_values(by=['TranscriptID', 'ExonRank']).reset_index(drop=True)
    
    # Add data to exontable
    for row in input_df.itertuples():
        exontable = pd.concat([
            exontable,
            pd.DataFrame({
                'GeneID': [row.GeneID],
                'TranscriptID': [row.TranscriptID],
                'ProteinID': [f"{row.TranscriptID}_PROTEIN"],
                'Strand': [row.Strand],
                'ExonID': [row.ExonID],
                'ExonRegionStart': [row.ExonRegionStart],
                'ExonRegionEnd': [row.ExonRegionEnd],
                'ExonRank': [row.ExonRank],
                'cDNA_CodingStart': [row.cDNA_CodingStart],
                'cDNA_CodingEnd': [row.cDNA_CodingEnd],
                'GenomicCodingStart': [row.GenomicCodingStart],
                'GenomicCodingEnd': [row.GenomicCodingEnd],
                'StartPhase': [row.StartPhase],
                'EndPhase': [row.EndPhase],
            })
        ], ignore_index=True)
    
    return exontable


def read_transcript_file(transcript_file):
    """
    Read the transcript CSV file and return a pandas DataFrame without modifications.
    """
    # Read data from CSV
    return pd.read_csv(transcript_file)


def add_to_tsl(input_df, tsl_df):
    """
    It adds the new transcripts to the TSL table downloaded from Ensembl.
    """
    subset = input_df[['Species', 'TranscriptID']].drop_duplicates()
    for row in subset.itertuples():
        tsl_df = tsl_df.append(
            {
                'Species':
                row.Species,
                'Name':
                row.TranscriptID,
                'TranscriptID':
                row.TranscriptID,
                'Source':
                'user',
                'ExperimentSource':
                'user',
                'Biotype':
                'protein_coding',
                'Flags':
                np.nan,
                'Version':
                str(datetime.datetime.now().isoformat(timespec='hours')),
            },
            ignore_index=True)
    return tsl_df


def add_sequences(input_df, seqrecords):
    """
    It modifies seqrecords by appending the new sequences.
    """
    for row in input_df.itertuples():
        name = f'{row.Species}:{row.GeneID}'
        seqrecords.append(
            SeqRecord(
                Seq(row.NucleotideSequence, SingleLetterAlphabet()),
                id=name,
                name=name,
                description=
                f'{row.Species}:{row.GeneID} {row.ExonID} na:na:na:{row.ExonRegionStart}:{row.ExonRegionEnd}:{row.Strand}' # pylint: disable=line-too-long
            ))

def main():
    """Main script function to add user transcript data."""
    args = parse_command_line().parse_args()
    paths = get_output_paths(args)

    print("Reading user input...")
    input_df = read_transcript_file(paths.input)  

    exontable = read_exon_file(paths.exontable)
    tsl_df = pd.read_csv(paths.tsl)

    print("Checking input...")
    check_input(input_df, exontable)

    print("Adding new transcripts...")
    new_exontable = add_to_exontable(input_df, exontable)

    new_tsl = add_to_tsl(input_df, tsl_df)

    seqrecords = read_fasta(paths.sequences)
    add_sequences(input_df, seqrecords)

    print("Saving...")
    _write_fasta(paths.sequences, seqrecords)
    new_tsl.to_csv(paths.tsl, index=False)
    new_exontable.to_csv(paths.exontable, sep='\t', index=False)

    print("Finished!")



if __name__ == '__main__':
    main()