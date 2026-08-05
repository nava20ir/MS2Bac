import pandas as pd
from io import StringIO
import numpy as np
"""
Usage: python prepare_from_casanovo.py <path-to-mztab-file> <Experiment_name>

"""
def read_mztab(mztab_file, skip_prefixes=("MTD",)):
    """
    Read an mzTab file into a pandas DataFrame while
    skipping lines that start with specified prefixes.

    Parameters
    ----------
    mztab_file : str
        Path to mzTab file.
    skip_prefixes : tuple
        Line prefixes to ignore.

    Returns
    -------
    pandas.DataFrame
    """
    with open(mztab_file) as f:
        lines = [
            line for line in f
            if not line.startswith(skip_prefixes)
        ]

    return pd.read_csv(
        StringIO("".join(lines)),
        sep="\t"
    )




if __name__ == '__main__':
    import sys
    df = read_mztab(sys.argv[1])
    df["score"] = df["opt_global_aa_scores"].apply(lambda x: np.mean([float(v) for v in x.split(",")]))
    df = (df.sort_values("score", ascending=False).drop_duplicates(subset="sequence", keep="first").reset_index(drop=True))
    df = df.rename(columns={"score":'Score', "sequence":"Sequence"})
    df = df[['Score','Sequence']]
    df['Experiment'] = sys.argv[2]
    df.to_csv('psm.csv',index=False)