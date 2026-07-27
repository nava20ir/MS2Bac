

from collections import defaultdict
import os
import numpy as np
import pandas as pd
from pathlib import Path


def combine_tsv_files(path, pattern="*.tsv"):
    """
    Combine all TSV files in a directory into a single DataFrame.

    Parameters
    ----------
    path : str or Path
        Directory containing TSV files.
    pattern : str, optional
        Glob pattern for file matching (default: '*.tsv').

    Returns
    -------
    pd.DataFrame
        Combined DataFrame.
    """
    path = Path(path)

    files = sorted(path.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No TSV files found in {path}")

    dfs = []
    for file in files:
        df = pd.read_csv(file, sep="\t")
        df["source_file"] = file.name  # optional: track origin
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)



def get_best_hit(df, by_col='hyperscore'):
    '''
    Get only the best peptide per psm

    Parameters
    ----------
    df : pandas DataFrame

    by_col: string , default hpyerscore
        column name on which best hit should be selected


    Returns
    --------
    # pandas Series:
        best hit for a scannumber based on by_col
    '''
    return df.sort_values(by=by_col, ascending=False).iloc[0]




def clean_msfragger_df(total):
    # Clean df
    total = total.rename(columns={'peptide': 'Sequence', 'hyperscore':'Score'})
    total['Sequence'] = total['Sequence'].str.replace('I|L', 'J', regex=True)
    total = total[['Experiment', 'Sequence', 'Score']]
    return total





def combine_msfragger_psms(path_input, experiment_name):
    '''
    Combine all MSfragger search outputs to one dataframe. Select only the best matching psm.

    Parameters
    ----------
    path_input : string
        path to MSfragger output files
    it: int
        first or second iteration


    Returns
    --------
    pd.Dataframe
        Combined MSfragger output file with only the best matching psm
        Writes additionally output files

    '''
    total = combine_tsv_files(path_input, pattern="*.tsv")
    total = total.groupby('scannum', as_index=False).agg(get_best_hit)
    total['Experiment'] = experiment_name
    total = clean_msfragger_df(total)
    print('Processing is done.')
    total.to_csv(os.path.join(path_input, 'psm.csv'), index=False)
        


if __name__ == "__main__":
    import sys
    path_input = sys.argv[1]  # Get the input path from command line arguments
    experiment_name = sys.argv[2]  # Get the experiment name from command line arguments
    combine_msfragger_psms(path_input, experiment_name)