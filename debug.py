"""This code serves for debugging the add_transcripts code.
Namely, we observed that some lines from the previous exonstable are removed after running this function.
This code serves to explain why does it happen and how to prevent it in the future.
"""

import pandas as pd

def filter_ensmus_rows(file1, file2, output1, output2):
    # Read the TSV files
    df1 = pd.read_csv(file1, sep='\t', header=None)
    df2 = pd.read_csv(file2, sep='\t', header=None)
    
    # Filter rows where the first column starts with 'ENSMUS'
    filtered_df1 = df1[df1[0].astype(str).str.startswith('ENSMUS')]
    filtered_df2 = df2[df2[0].astype(str).str.startswith('ENSMUS')]
    
    # Save the filtered data to new TSV files
    filtered_df1.to_csv(output1, sep='\t', index=False, header=False)
    filtered_df2.to_csv(output2, sep='\t', index=False, header=False)
    
    print(f'Filtered files saved as {output1} and {output2}')

def cast_floats_to_ints(file, output):
    df = pd.read_csv(file, sep='\t', header=None)
    df = df.applymap(lambda x: int(x) if isinstance(x, float) else x)
    df.to_csv(output, sep='\t', index=False, header=False)
    print(f'Floats casted to integers in {output}')

def compare_files(file1, file2, same_output, diff_output):
    df1 = pd.read_csv(file1, sep='\t', header=None, dtype=str)
    df2 = pd.read_csv(file2, sep='\t', header=None, dtype=str)
    
    common_rows = df1.merge(df2, how='inner')
    different_rows = pd.concat([df1, df2]).drop_duplicates(keep=False)
    
    common_rows.to_csv(same_output, sep='\t', index=False, header=False)
    different_rows.to_csv(diff_output, sep='\t', index=False, header=False)
    
    print(f'Same lines saved to {same_output}')
    print(f'Different lines saved to {diff_output}')

def check_empty_values(file):
    df = pd.read_csv(file, sep='\t', header=None, dtype=str)
    empty_columns = df.isin(['']).sum()
    empty_columns = empty_columns[empty_columns > 0]
    if empty_columns.empty:
        print(f'No empty values found in {file}')
    else:
        print(f'Empty values found in columns of {file}:')
        print(empty_columns.to_string())

# Example usage
file1 = '/Users/milicaurosevic/TP53_3/exonstable_copy.tsv'  # before add_transc 
file2 = '/Users/milicaurosevic/TP53/Ensembl/exonstable.tsv'  # Replace with your actual file name
output1 = 'MouseBefore.tsv'
output2 = 'MouseAfter.tsv'

filter_ensmus_rows(file1, file2, output1, output2)
cast_floats_to_ints(output2, 'MouseAfterCasted.tsv')
compare_files('MouseBefore.tsv', 'MouseAfterCasted.tsv', 'SameLines.tsv', 'DifferentLines.tsv')
check_empty_values('SameLines.tsv')
check_empty_values('DifferentLines.tsv')