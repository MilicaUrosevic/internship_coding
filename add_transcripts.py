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
    'ProteinID',
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
    # Normalize input table column names
    input_columns = input_df.columns.str.strip().str.lower()
    
    # Normalize expected column names
    expected_columns = {col.lower() for col in COLUMN_NAMES}
    
    print(f"Normalized input columns: {list(input_columns)}")
    print(f"Normalized expected columns: {list(expected_columns)}")
    
    # Check for missing columns
    for col in expected_columns:
        if col not in input_columns:
            raise Exception(f'The {col} column is missing in the input table.')
    
    # Warn about unused columns
    for col in input_columns:
        if col not in expected_columns:
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



#  this is fine, only checks exon id
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
    for index, row in input_df.iterrows():  # Use iterrows to loop through rows with their index
        if row.ExonRegionStart > row.ExonRegionEnd:
            raise Exception(
                f"Row {index}: ExonRegionStart ({row.ExonRegionStart}) should be less than ExonRegionEnd ({row.ExonRegionEnd})."
            )
        if row.GenomicCodingStart >= row.GenomicCodingEnd:
            raise Exception(
                f"Row {index}: GenomicCodingStart ({row.GenomicCodingStart}) should be less than GenomicCodingEnd ({row.GenomicCodingEnd})."
            )
        if len(row.NucleotideSequence) != (row.ExonRegionEnd - row.ExonRegionStart + 1):
            raise Exception(
                f"Row {index}: The sequence length estimated using ExonRegionStart ({row.ExonRegionStart}) "
                f"and ExonRegionEnd ({row.ExonRegionEnd}) disagrees with the actual NucleotideSequence length "
                f"({len(row.NucleotideSequence)} vs {row.ExonRegionEnd - row.ExonRegionStart + 1})."
            )


# this function checks only the csv table
def check_input(input_df, exontable):

  """  It raises an error if there is a problem with the user-defined transcripts. """
  _check_col_names(input_df)
  _check_species_name(input_df)
  _check_transcript(input_df, exontable)
  _check_pair_order(input_df)
  _check_exon(input_df, exontable)


def track_changes(file_path, description, delimiter=',', sep='\t'):
    """Check and print the number of rows in a given file before and after operations."""
    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path, delimiter=delimiter)
        elif file_path.endswith('.tsv'):
            df = pd.read_csv(file_path, sep=sep)
        else:
            print(f"[WARNING] Unknown format for {file_path}, skipping tracking.")
            return
        
        print(f"{description}: {file_path} -> {len(df)} rows")
    except Exception as e:
        print(f"[ERROR] Could not read {file_path}: {e}")



#  THIS IS SUSPICIOUS, actually if GenomicCodingEnd or GenomicCodingEnd is NaN, everything else is supposed to be also NaN that includes them
"""def add_to_exontable(input_df, exontable):
    """ """
    It adds the new exons and transcripts to the Ensembl's exontable.
    also check these calculations
    """ """
    transcript_id = ""
    cdna_start = 0
    cdna_end = 0
    for row in input_df.itertuples():
        if transcript_id != row.TranscriptID or cdna_end == 0:
            transcript_id = row.TranscriptID
            coding_len = 0 if pd.isna(row.GenomicCodingStart) or pd.isna(row.GenomicCodingEnd) else row.GenomicCodingEnd - row.GenomicCodingStart
            exon_len = row.ExonRegionEnd - row.ExonRegionStart
            cdna_start = exon_len - coding_len + 1 if coding_len > 0 else np.nan
            cdna_end = cdna_start + coding_len if coding_len > 0 else np.nan
        else:
            cdna_start = cdna_end + 1
            cdna_end += (row.GenomicCodingEnd - row.GenomicCodingStart + 1)
        exontable = pd.concat([exontable, pd.DataFrame([{
            'GeneID': row.GeneID,
            'TranscriptID': row.TranscriptID,
            'ProteinID': row.ProteinID if pd.notna(row.ProteinID) else np.nan,
            'Strand': row.Strand,
            'ExonID': row.ExonID,
            'ExonRegionStart': row.ExonRegionStart,
            'ExonRegionEnd': row.ExonRegionEnd,
            'ExonRank': row.ExonRank,
            'cDNA_CodingStart': cdna_start,
            'cDNA_CodingEnd': cdna_end,
            'GenomicCodingStart': row.GenomicCodingStart if pd.notna(row.GenomicCodingStart) else np.nan,
            'GenomicCodingEnd': row.GenomicCodingEnd if pd.notna(row.GenomicCodingEnd) else np.nan,
            'StartPhase': row.StartPhase,
            'EndPhase': row.EndPhase,
        }])], ignore_index=True)

    return exontable"""

def add_to_exontable(input_df, exontable):
    print(f"Before adding: {len(exontable)} rows")  # Broj redova pre dodavanja
    
    transcript_id = ""
    cdna_start = 0
    cdna_end = 0
    for row in input_df.itertuples():
        if transcript_id != row.TranscriptID or cdna_end == 0:
            transcript_id = row.TranscriptID
            coding_len = 0 if pd.isna(row.GenomicCodingStart) or pd.isna(row.GenomicCodingEnd) else row.GenomicCodingEnd - row.GenomicCodingStart
            exon_len = row.ExonRegionEnd - row.ExonRegionStart
            cdna_start = exon_len - coding_len + 1 if coding_len > 0 else np.nan
            cdna_end = cdna_start + coding_len if coding_len > 0 else np.nan
        else:
            cdna_start = cdna_end + 1
            cdna_end += (row.GenomicCodingEnd - row.GenomicCodingStart + 1)
        
        if pd.isna(row.ExonID) or pd.isna(row.TranscriptID):
            print(f"[WARNING] Skipping row because ExonID or TranscriptID is NaN: {row}")

        
        before_concat = len(exontable)
        exontable = pd.concat([exontable, pd.DataFrame([{
            'GeneID': row.GeneID,
            'TranscriptID': row.TranscriptID,
            'ProteinID': row.ProteinID if pd.notna(row.ProteinID) else np.nan,
            'Strand': row.Strand,
            'ExonID': row.ExonID,
            'ExonRegionStart': row.ExonRegionStart,
            'ExonRegionEnd': row.ExonRegionEnd,
            'ExonRank': row.ExonRank,
            'cDNA_CodingStart': cdna_start,
            'cDNA_CodingEnd': cdna_end,
            'GenomicCodingStart': row.GenomicCodingStart if pd.notna(row.GenomicCodingStart) else np.nan,
            'GenomicCodingEnd': row.GenomicCodingEnd if pd.notna(row.GenomicCodingEnd) else np.nan,
            'StartPhase': row.StartPhase,
            'EndPhase': row.EndPhase,
    }])], ignore_index=True)

    after_concat = len(exontable)
    if after_concat == before_concat:
        print(f"[WARNING] Row with ExonID={row.ExonID} and TranscriptID={row.TranscriptID} was not added!")


    print(f"After adding: {len(exontable)} rows")  # Broj redova posle dodavanja
    print(f"Before returning from add_to_exontable: {len(exontable)} rows")

    return exontable


# CHECK THIS ALSO
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
    
    print("Checking initial file sizes...")
    track_changes(paths.input, "Initial input file")
    track_changes(paths.exontable, "Initial exon table")
    track_changes(paths.tsl, "Initial TSL table")
    track_changes(paths.sequences, "Initial FASTA sequences")
    
    print("Reading user input...")
    input_df = pd.read_csv(paths.input, delimiter=',')  # changed
    # Proveravamo koliko redova zapravo ima u fajlu

    # Proveravamo koliko redova stvarno ima u fajlu pre nego što ga pročita read_exon_file
    df_check = pd.read_csv(paths.exontable, sep='\t')
    print(f"Directly reading exonstable.tsv (before read_exon_file): {len(df_check)} rows")

    # Sada koristimo read_exon_file i proveravamo da li on skraćuje tabelu
    exontable = read_exon_file(paths.exontable)
    print(f"After read_exon_file: {len(exontable)} rows")

    tsl_df = pd.read_csv(paths.tsl)
    seqrecords = read_fasta(paths.sequences)

    track_changes(paths.input, "After reading input file")
    track_changes(paths.exontable, "After reading exon table")
    track_changes(paths.tsl, "After reading TSL table")
    track_changes(paths.sequences, "After reading FASTA sequences")

    print("Checking input...")
    check_input(input_df, exontable)
    track_changes(paths.input, "After check_input")

    print("Adding new transcripts...")
    new_exontable = add_to_exontable(input_df, exontable)
    track_changes(paths.exontable, "After add_to_exontable")

    new_tsl = add_to_tsl(input_df, tsl_df)
    track_changes(paths.tsl, "After add_to_tsl")

    print("Modifying sequences...")
    add_sequences(input_df, seqrecords)
    track_changes(paths.sequences, "After add_sequences")

    print("Saving...")
    _write_fasta(paths.sequences, seqrecords)
    
    print(f"Before saving exonstable.tsv: {len(new_exontable)} rows")
    print("NaN counts per column before saving:")
    print(new_exontable.isnull().sum())  # Proveri koliko ima NaN vrednosti u svakoj koloni
    print("Checking for duplicate rows before saving:")
    print(new_exontable.duplicated().sum())  # Proveri broj duplikata

    print("NaN counts per column before saving:")
    print(new_exontable.isnull().sum())  # Proveri da li postoje NaN vrednosti u nekoj koloni


    new_tsl.to_csv(paths.tsl, index=False)
    new_exontable.to_csv(paths.exontable, sep='\t', index=False, na_rep="NA")

    track_changes(paths.tsl, "After saving TSL table")
    track_changes(paths.exontable, "After saving exon table")
    track_changes(paths.sequences, "After saving FASTA sequences")

    print("Finished!")


if __name__ == '__main__':
    main()
