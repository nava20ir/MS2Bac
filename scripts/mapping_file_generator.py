#*******************************************************************************
# Copyright 2024 Miriam Abele <m.abele@tum.de> and Christina Ludwig (tina.ludwig@tum.de)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#*******************************************************************************

####################
##### Packages #####
####################

import socket
import sys
from collections import defaultdict
import os
import tqdm

from collections import defaultdict
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

import pyteomics.fasta as fasta



#####################
##### Functions #####
#####################

def create_iDB(path_digest, path_taxon_annot, it):
    '''
    Function builts up a bacci database required for identification.
    This process takes some time, but can be re-used any time.
    I and L will be treated equally, because they cannot be distinguised by
    MS due to their equal mass.


    Parameters
    ----------
    path_digest: string
        path to folder with digested files
    path_taxon_annot : string
        ath to excel file with all NCBI taxon info

    Returns
    --------
    idb : defaultdict
        identification database;
        sequences as keys and taxonids as values, can afterwards be searched.
    bacci_finder : dict
        translate taxon_id in organism name (o_name)

    '''

    # Variable assignments
    idb = defaultdict(set)
    digests = os.listdir(path_digest)

    if digests == 0:
        print('############## Be aware that your digest folder is empty! ############')

    bacci_finder = dict()
    annotation = pd.read_csv(path_taxon_annot)

    # create defaultdict with sequences as keys and taxonids as values, taxonid vary depending on input
    for bacci in tqdm.tqdm(digests, mininterval=10):
        try:
            colname = get_accession_assembly_column(annotation)
            bacci_seqs = set(np.genfromtxt(os.path.join(path_digest, bacci, 'protein_seq.csv'), dtype='str', skip_header=1))
            len_digest = len(bacci_seqs)

            ncbi_id = bacci.split('_')[1]
            genus_tax_id = annotation[annotation[colname] == ('GCF_'+ncbi_id)]['genus'].item()

            if it == 1:
                taxon_id = annotation[annotation[colname] == ('GCF_'+ncbi_id)]['species'].item()
                o_name = annotation[annotation[colname] == ('GCF_'+ncbi_id)]['organism_name'].item()
            elif it == 2:
                taxon_id = annotation[annotation[colname] == ('GCF_'+ncbi_id)][colname].item()
                annotation['o_name_secit'] = annotation['organism_name'] + '; ' + annotation['infraspecific_name']
                o_name = annotation[annotation[colname] == ('GCF_'+ncbi_id)]['o_name_secit'].item()
            bacci_finder[taxon_id] = (o_name, genus_tax_id, len_digest)

        except ValueError:
            print('No annotation available. This should not be the case, because all of them were downloaded! Please check!')
            continue
        except OSError:
            print('OS Error occured, probably due to another folder. Continue')
            continue


    # iterate through sequences of bacci from bacci
        for seq in bacci_seqs:
            seq = seq.replace('I', 'J').replace('L', 'J')
            idb[seq].add(taxon_id)


    print(f'Length of the idb is {len(idb)}')
    return idb, bacci_finder


def get_accession_assembly_column(df):
    '''
        Function returns the column name of the assembly accession column in the assembly_summary.txt file

        Parameters
        ----------
        df : pd.DataFrame
            assembly_summary.txt file from NCBI

        Returns
        --------
        string
            column name of the assembly accession column
    '''

    for col in df.columns:
        if 'assembly_accession' in col:
            return col

if __name__ == "__main__":
    import sys
    root = sys.argv[1]
    path_digest = os.path.join(root, '03_digest')
    path_taxon_annot = os.path.join(root, '00_metafiles', 'download_identification_database.csv')
    idb, bacci_finder = create_iDB(path_digest, path_taxon_annot, 1)

    import pickle
    import gzip

    with gzip.open("idb.pkl.gz", "wb", compresslevel=1) as f:
        pickle.dump(idb, f, protocol=pickle.HIGHEST_PROTOCOL)

    with gzip.open("bacci_finder.pkl.gz", "wb", compresslevel=1) as f:
        pickle.dump(bacci_finder, f, protocol=pickle.HIGHEST_PROTOCOL)

        
# how to open the database:
# with gzip.open("idb.pkl.gz", "rb") as f: idb = pickle.load(f)
# with gzip.open("bacci_finder.pkl.gz", "rb") as f: bacci_finder = pickle.load(f)