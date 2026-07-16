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

# The script downloads proteomic information from NCBI (citation: National Center for Biotechnology Information (NCBI)[Internet]. Bethesda (MD): National Library of Medicine (US),
# National Center for Biotechnology Information; [1988] - [cited 2024 Jan 23]. Available from https://www.ncbi.nlm.nih.gov/).
# Please respect the terms and conditions from NCBI!


# system packages/folder structure/download
import os
import tqdm
import argparse
import wget
import gzip
from datetime import datetime, date
import time
import shutil
import urllib
import urllib.request as request
from contextlib import closing

# data manipulation packages
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from collections import defaultdict

# proteomics packages
import pyteomics
from pyteomics import fasta
from pyteomics.parser import cleave, expasy_rules
from ete3 import NCBITaxa


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

######################## functions for metadata pre-processsing ####################################
  
def taxon_annotation(assembly, taxon, ncbi):
    ''' 
        Function creates annotation dictionaries

        Parameters
        ----------
        assembly : string
            path to the assembly_summary.txt file from NCBI.
        taxon : string
            taxon of interest, e.g. 'phylum' or 'genus'

        ncbi: NCBITaxa object from the ete3 package
            
        Returns
        --------
        taxon_dict: dictionary
            Mapper for NCBI taxid (key) to NCBI taxon identifier (value)
        name_dict: dictionary
            Mapper for NCBI taxid (key) to NCBI taxon name (value)
    '''

    taxon_dict = dict()
    name_dict = dict()

    for t_id in tqdm.tqdm_notebook(assembly.taxid):
        try:
            lineage = ncbi.get_lineage(int(t_id))
            rank = ncbi.get_rank(lineage)

            for key, value in rank.items():
                if value == taxon:
                    name =  ncbi.get_taxid_translator([key])[key]
                    taxon_dict[t_id] = key
                    name_dict[t_id] = name
                    break
                else:
                    continue

        except IndexError:
            print(ncbi.get_lineage(taxon))
    
        except ValueError:
            # For some taxons 
            taxon_dict[t_id] = np.nan
            name_dict[t_id] = np.nan 


    return taxon_dict, name_dict


def metadata_preprocessing(assembly_path, ani_path):
    '''
        Function cleans assembly_summary.txt file and adds annotation. Additionally it merges the ANI from the ANI_report.txt file. 

        Parameters
        ----------
        assembly_path : string
            path to the assembly_summary.txt file from NCBI.
        proteome_path : string
            path to the ANI_report_prokaryotes.txt file


        Returns
        --------
        no return, but saves a new file named xxx to 00_metafiles
    '''
    # Load file
    assembly_summary = pd.read_csv(assembly_path, sep='\t', low_memory=False, skiprows=1,error_bad_lines=False, encoding='latin-1')

    # Clean file
    assembly_summary['ftp_path'] = assembly_summary['ftp_path'].replace('na', np.nan)
    assembly_summary = assembly_summary.dropna(subset = ['ftp_path'])

    # Update ncbi taxon database via ncbi FTP server
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()

    # Taxon Annotation
    taxon_dict = dict()
    name_dict = dict()
    for taxon in ['phylum', 'genus', 'species']:
        taxon_dict, name_dict = taxon_annotation(assembly= assembly_summary, 
                                             taxon=taxon, ncbi=ncbi)
        assembly_summary[taxon] = assembly_summary['taxid'].map(taxon_dict)
        assembly_summary['ncbi_'+str(taxon)+'_name']= assembly_summary['taxid'].map(name_dict)

    # Save file
    output_path = assembly_path.replace('.txt', '_annot.csv')
    assembly_summary.to_csv(output_path)

    # Print summary
    print(assembly_summary[['phylum', 'genus', 'species']].nunique())

    # Load ANI file
    ani = pd.read_csv(ani_path,sep='\t')  

    # Map to assembly_summary file via refseq identifier
    df_all = assembly_summary.merge(ani, left_on=get_accession_assembly_column(assembly_summary), right_on ='refseq-accession', how='left')

    output_path = assembly_path.replace('.txt', '_annot_withANI.csv')
    df_all.to_csv(output_path)

####################### functions for metadata pre-processsing end ###################################

def ncbi_download_reference(
    assembly_summary_path,
    proteome_path,
    metadata_output_path,
    timeout=60,
    max_attempts=3
):
    """
    Download one representative/reference proteome per species from NCBI.
    """

    # ------------------------------------------------------------------
    # Load metadata
    # ------------------------------------------------------------------
    assembly_summary = pd.read_csv(
        assembly_summary_path,
        low_memory=False
    )

    print(
        f"Number of proteomes before filtering: "
        f"{len(assembly_summary)}"
    )

    assembly_summary = assembly_summary[
        assembly_summary["taxonomy-check-status"] == "OK"
    ]

    assembly_summary = assembly_summary[
        assembly_summary["refseq_category"].isin(
            ["representative genome", "reference genome"]
        )
    ]

    print(
        f"Number of proteomes after filtering: "
        f"{len(assembly_summary)}"
    )

    accession_col = get_accession_assembly_column(
        assembly_summary
    )

    reference = assembly_summary[
        ["ftp_path", "species", accession_col]
    ]

    print(
        f"Current database contains {len(reference)} "
        f"reference/representative genomes"
    )

    # ------------------------------------------------------------------
    # Download
    # ------------------------------------------------------------------
    existing_files = set(os.listdir(proteome_path))

    species_download = {}
    failed_urls = []

    for download_path, species_id, accession in tqdm.tqdm(
        reference.values,
        mininterval=10
    ):

        # keep only one genome per species
        if species_id in species_download:
            continue

        species_download[species_id] = accession

        #organism_id = download_path.split("/")[-1]
        organism_id = os.path.basename(download_path.rstrip("/"))
        outfile = os.path.join(
            proteome_path,
            organism_id + "_protein.faa.gz"
        )

        if os.path.basename(outfile) in existing_files:
            continue

        url = (
            f"{download_path}/"
            f"{organism_id}_protein.faa.gz"
        )

        print("\n" + "=" * 80)
        print(f"Species   : {species_id}")
        print(f"Accession : {accession}")
        print(f"URL       : {url}")
        print("=" * 80)

        success = False

        for attempt in range(1, max_attempts + 1):

            try:

                print(
                    f"Attempt {attempt}/{max_attempts}"
                )

                req = request.Request(
                    url,
                    headers={
                        "User-Agent":
                        "MS2Bac/1.0 (NCBI protein download)"
                    }
                )

                with closing(
                    request.urlopen(
                        req,
                        timeout=timeout
                    )
                ) as response:

                    with open(outfile, "wb") as f:
                        shutil.copyfileobj(
                            response,
                            f
                        )

                print(
                    f"SUCCESS: {organism_id}"
                )

                success = True
                break

            except urllib.error.HTTPError as err:

                print("\nHTTP ERROR")
                print(
                    f"Code   : {err.code}"
                )
                print(
                    f"Reason : {err.reason}"
                )
                print(
                    f"URL    : {url}"
                )

                # File not present on NCBI
                if err.code == 404:
                    print(
                        "Protein file not found. "
                        "Skipping."
                    )
                    break

            except urllib.error.URLError as err:

                print("\nURL ERROR")
                print(
                    f"Reason : {err.reason}"
                )
                print(
                    f"URL    : {url}"
                )

            except Exception as err:

                print("\nUNEXPECTED ERROR")
                print(
                    f"Type   : {type(err).__name__}"
                )
                print(
                    f"Error  : {err}"
                )
                print(
                    f"URL    : {url}"
                )

            # Wait before retry
            time.sleep(2)

        if not success:
            failed_urls.append(url)

    # ------------------------------------------------------------------
    # Save metadata
    # ------------------------------------------------------------------
    selected_accessions = set(
        species_download.values()
    )

    ncbi_species_df = assembly_summary[
        assembly_summary[accession_col].isin(
            selected_accessions
        )
    ]

    output_file = os.path.join(
        metadata_output_path,
        "download_identification_database.csv"
    )

    ncbi_species_df.to_csv(
        output_file,
        index=False
    )

    print(
        f"\n{len(ncbi_species_df)} proteomes "
        f"written to metadata file"
    )

    print(
        f"{len(failed_urls)} downloads failed"
    )

    if failed_urls:
        print("\nFailed URLs:")
        for url in failed_urls[:20]:
            print(url)

    return failed_urls
    

########################## functions for ncbi download #############################################

def ncbi_download_referenceOLD(assembly_summary_path, proteome_path, metadata_output_path):
    '''
        Function downloads bacteria proteomes from NCBI.
        Filtering critera are: taxonomy == ok, assembly_level == 'representative' or 'reference

        Parameters
        ----------
        assembly_summary_path : string
            path to assembly summary file from NCBI.
        proteome_path : string
            path to the folder location where proteomes are stored
        metadata_output_path : string
            path to the folder location where download metadata should be stored

        Returns
        --------
        no return
    '''

    # load assebmly summary file that includes annotation and ani results from NCBI
    assembly_summary = pd.read_csv(assembly_summary_path, low_memory=False)

    # filter for taxonomy == 'ok' and assembly_level == 'representative' or 'reference
    print(f'Number of proteomes in assembly summary before filtering: {len(assembly_summary)}\n')

    assembly_summary = assembly_summary[assembly_summary['taxonomy-check-status'] == 'OK']
    assembly_summary = assembly_summary[(assembly_summary['refseq_category'] == 'representative genome') | (assembly_summary['refseq_category'] == 'reference genome')]

    print(f'Number of proteomes in assembly summary after filtering: {len(assembly_summary)}\n')


    print(f'Number of proteomes in assembly summary after filtering sp.: {len(assembly_summary)}\n')

    print(f'Current ncbi database contains {len(assembly_summary)} complete genomes')


    # download from ncbi, download requires a 2 sec sleep in order to continue
    species_download = {}
    reference = assembly_summary[['ftp_path', 'species', get_accession_assembly_column(assembly_summary)]]

    for download_path, ncbi_taxon_id, accession in tqdm.tqdm(reference.values, mininterval=10):

        # Check if the species already exists:
        if ncbi_taxon_id not in species_download:
            species_download[ncbi_taxon_id] = (accession)

            # Check if file already exists
            organism_id = download_path.split('/')[-1]

            if (organism_id+'_protein.faa.gz') in os.listdir(proteome_path):
                continue

            else:
                # download from NCBI
                u = (download_path+'/'+organism_id+'_protein.faa.gz')
                print(u)
                time.sleep(2)

                try:
                    attempts = 0

                    while attempts < 3:
                        print(f'{attempts+1}\n')

                        try:
                            print(f'Start request {attempts+1}')
                            with closing(request.urlopen(u, timeout=2)) as r:
                                print('Open file')
                                with open((proteome_path+organism_id+'_protein.faa.gz'), 'wb') as f:
                                    print('move file')
                                    shutil.copyfileobj(r, f)
                                    print('Done. Next File!')
                                    break

                        except (urllib.error.HTTPError, urllib.error.URLError) as error:
                            attempts += 1
                            print(type(error))


                    print('End download')

                except ValueError:
                    print('Exception occured! Please try again\n\n')
                    pass

        else:
            pass

    # print new metadata_species_only to csv
    ncbi_species_df = assembly_summary[assembly_summary[get_accession_assembly_column(assembly_summary)].apply(lambda i: i in species_download.values())]
    print(f'{len(ncbi_species_df)} proteomes in metadata_species file!')
    ncbi_species_df.to_csv(os.path.join(metadata_output_path, 'download_identification_database.csv'))


########################## functions for ncbi download end #############################################



########################## functions for unzipping files #############################################

# Unzip files and store them in a folder called unzip
def unzip_files(proteome_path, unzip_path):
    '''
    Proteomes from NCBI download are stored as .gz files. This function unzips them.

    Parameters
    ----------
    proteome_path : string
        path to the folder location where proteomes are stored
    unzip_path : string
        path to the folder location where unzipped proteomes are stored


    Returns
    --------
    output is stored in unzip_path specified in input

    '''

    counter=0
    proteome_bacci = [s for s in (os.listdir(proteome_path))]

    print('Start unzipping files!')
    for i, file in enumerate(proteome_bacci):
        counter+=1
	
        if (file + '.fasta') in os.listdir(unzip_path):
            continue

        else:
            if 'gz' in file:
                try:
                    input = gzip.GzipFile(proteome_path + '/' + file, 'rb')
                    s = input.read()
                    input.close()

                    output = open(unzip_path + '/' + file + ".fasta", 'wb')
                    output.write(s)
                    output.close()
                except EOFError:
                    print(file)
            else:
                continue
                
    print('All files are unzipped and ready for digest!')


########################## functions for digest #############################################

def protein_factory(rule, min_length, missed_cleavages):
    '''
    In silico digest of proteins with pyteomics package

    Parameters
    ----------
    rule : string
        enzyme for digest
    min_length : int
        filter for minimal peptide length
    missed_cleavages: into
        number of allowed missed cleavages

    Returns
    --------
    set
        product peptides
    '''

    r = expasy_rules[rule]
    return lambda protein: cleave(protein, r, missed_cleavages, min_length)


def split_description(string):
    '''
    plit string according to . and whitespace character and returns first element

    Parameters
    ----------
    string : string
        protein

    Returns
    --------
    string
        first element of the split
    '''

    string = string.replace('.', ' ')
    return string.split(" ")[1]


def digest_peptides_df(path, missed_cleavages=2, rule="trypsin", min_length=7, with_protein=False):
    '''
    Performs an in silico digests of proteomes with filtering criteria and stores them in a peptide_seq dataframe

    Parameters
    ----------
    path : string
        path to one specific proteome in the unzip folder
    mc : int
        missed_cleavages, default 2
    rule : string
        enzyme for digest, default "trypsin"
    min_length : int
        filter for minimal peptide length, default 7
    with_protein : bool
        Add protein name to digest file, default False

    Returns
    --------
    pd.DataFrame
        column identifier and column peptide_seq
    '''
    digest = protein_factory(rule, min_length, missed_cleavages)

    with open(path, 'rt') as f:
        fasta = pyteomics.fasta.read(f)

        protein_a = []
        peptide_a = []
        for protein in fasta:
            peptides = list(digest(protein.sequence))
            peptides = [
                p
                for p in peptides
                if (
                    (len(p) <= 30)
                    & ("X" not in p)
                    & ("U" not in p)
                    & ("B" not in p)
                    & ("i" not in p)
                    & ("Z" not in p)
                )
            ]
            protein_a.extend(
                [split_description(protein.description) for _ in range(len(peptides))]
            )
            peptide_a.extend(peptides)

        return pd.DataFrame({"identifier": protein_a, "peptide_seq": peptide_a})


def write_f(output_folder, df):
    '''
    Function writes digest unzip_files

    Parameters
    ----------
    output_folder : string
        file location where digested files should be stored


    Returns
    --------
    writes protein_seq.csv and protein_all.csv files to output_folder
    '''


    df["peptide_seq"].to_csv(
        os.path.join(output_folder, "protein_seq.csv"), index=False
    )
    df.to_csv(os.path.join(output_folder, "protein_all.csv"), index=False)


def run_fasta(path, output_folder, mc):
    '''
    Function writes digest unzip_files

    Parameters
    ----------
    path : string
        folder location to one specific proteome in the unzipped folder
    output_folder: string
        older location to one specific digest in the digest folder
    mc : int
        number of missed cleavages


    Returns
    --------
    writes files to output_folder
    '''

    df = digest_peptides_df(path, missed_cleavages=mc)
    os.mkdir(output_folder)
    write_f(output_folder, df)

def run_digest(mc, unzip_path, digest_path):
    '''
    Function digests all proteomes that are located in the unzip_path to peptides

    Parameters
    ----------
    unzip_path : string
        path to folder location where unzipped files are stored
    digest_path: string
        path to folder location where digests should be stored
    mc : int
        number of missed cleavages


    Returns
    --------
    writes digested files to digest_path
    '''

 
    print('Start digest!')
    for i in os.listdir(unzip_path):
        if (str(i)) in os.listdir(digest_path):
            continue
        else:
            try:
                if ('fasta' in i) and ('DNA' not in i):
                    run_fasta(unzip_path+'/'+str(i), (digest_path+'/'+str(i)+'/'), mc)
            except FileExistsError:
                print(i)
                continue
    print('All in silico digests performed!')


########################## functions for digest end #############################################




########################## functions for writing fasta #############################################
def write_fasta_2(digest_path, fasta_path, n_fastas = 4):
    '''
        Function writes two fasta files, splits all peptides evenly
        Every peptide is one 'fake' protein

        Parameters
        ----------
        digest_path: string
            path to folder location where digests should be stored
        fasta_path : string
            path to folder location where fastas should be store


        Returns
        --------
        writes fasta files to fasta_path
    '''

    files = os.listdir(digest_path)

    db = defaultdict(list)
    organism_id_dict = dict()

    # iterate through all reference proteome files that were digested
    print('Last step: Creation of fasta for search engine input starts!')
    for i, file_name in tqdm.tqdm(enumerate(files), mininterval=10):
        try:
            organism_seqs = np.genfromtxt(os.path.join(digest_path,file_name,'protein_seq.csv'), dtype='str', skip_header=1)
        except OSError:
            print(f'OS Error file: {file_name}')
            continue

        # We don't want to have duplicates --> set
        organism_set = set(organism_seqs)

        # save uniprot_id in organism dict and translate into a int value to reduce memory
        try:
            nacbi_id = file_name.split('_')[0]+'_'+file_name.split('_')[1]+'_'+file_name.split('_')[2]
            organism_id_dict[i] = nacbi_id
        except IndexError:
            nacbi_id

        # iterate through sequences and append db
        for seq in organism_set:
            db[seq].append(i)

    # Distribute all sequences equally to n fasta files
    fasta_entries = {i:list() for i in range(n_fastas)}
    for i, key in tqdm.tqdm(enumerate(db.keys()), mininterval=10):
        fasta_entries[i % n_fastas].append((str(i), key))

    # Write fastas
    for i, seq_lst in fasta_entries.items():
        date = str(datetime.now())[:10]
        file_name = f'{date}_{i}.fasta'
        path = os.path.join(fasta_path, file_name)
        fasta.write(seq_lst, path)
        print(f'Fasta {i} done!')

    print('Fastas are written. You can start now searching your MS data with those files!')
    print('Good luck with your identification :) ')

    del fasta_entries
    del db

######################### functions for writing fasta end ###########################################



###################### functions for second iteration download ######################################

def create_all_files(metadata_path, proteome_path, unzip_path, digest_path, fasta_path):
    '''
    function creates all folder required for the second iteration download

    Parameters
    ----------
    metadata_path : string
        path to metadata for second iteration
    proteome_path : string
        path to downloaded proteomes for second iteration
    unzip_path : string
        path to unzipped files for second iteration
    digest_path : string
        path to digest for second iteration
    fasta_path : string
        path to fasta for second iteration

    Returns
    --------
        makes all the directories.
    '''

    try:
        os.mkdir(metadata_path)
        os.mkdir(proteome_path)
        os.mkdir(unzip_path)
        os.mkdir(digest_path)
        os.mkdir(fasta_path)
        print('Genus-specific background database directories made. ')

    except OSError:
        ('Creation of directy failed. Maybe the folder already exists')


def ncbi_download_secit(root, goi, root_taxon_specific, p_id, proteome_path, digest_path):
    '''
    Function downloads bacteria proteomes from NCBI. It will only download the first strain per species which is in the metadatafile!

    Parameters
    ----------
    root : string
        path to metadata from ncbi
    goi : string
        genus of interest
    root_taxon_specific : string
        path to genus specific file location
    p_id : string
        project identifier

    Returns
    --------
        makes all the directories.
    '''

    # load metadata file
    root_metadata = os.path.dirname(root)
    metadata_path = os.path.dirname(root) + '/00_metafiles/assembly_summary_annot_withANI.csv'
    metadata = pd.read_csv(metadata_path, low_memory=False)

    # Filter for genus ncbi ID
    metadata = metadata[metadata['genus'] == int(goi)]
    metadata = metadata[metadata['taxonomy-check-status'] == 'OK']

    reference = metadata[['ftp_path', 'organism_name', 'infraspecific_name', get_accession_assembly_column(metadata)]]


    species_download= {}
    print(f'Current ncbi database contains {len(metadata)} genomes of the NCBI genus {goi}.')

    # download from ncbi, download requires a 8 sec sleep in order to continue

    for download_path, ncbi_taxon_id, strain_name, accession in tqdm.tqdm(reference.values, mininterval=10):
        # Check if the species already exists:
        real_name = ncbi_taxon_id + str(strain_name)

        if real_name not in species_download:
            species_download[real_name] = (accession)
            # Check if file already exists
            organism_id = download_path.split('/')[-1]

            if (organism_id+'_protein.faa.gz') in os.listdir(proteome_path):
                continue

            else:
                u = (download_path+'/'+organism_id+'_protein.faa.gz')
                print(u)
                time.sleep(2)

                try:
                    attempts = 0

                    while attempts < 3:
                        print(f'{attempts+1}\n')

                        try:
                            print(f'Start request {attempts+1}')
                            with closing(request.urlopen(u, timeout=2)) as r:
                                print('Open file')
                                with open((proteome_path+'/'+organism_id+'_protein.faa.gz'), 'wb') as f:
                                    print('move file')
                                    shutil.copyfileobj(r, f)
                                    print('Done. Next File!')
                                    break

                        except (urllib.error.HTTPError, urllib.error.URLError) as error:
                            attempts += 1
                            print(type(error))


                    print('End download')

                except ValueError:
                    print('Exception occured! Please try again\n\n')
                    pass

                except timeout:
                    shutil.rmtree(root_taxon_specific + '/00_metafiles/' + goi)
                    shutil.rmtree(root_taxon_specific + '/01_proteomes/' + goi)
                    shutil.rmtree(root_taxon_specific + '/02_unzipped_files/' + goi)
                    shutil.rmtree(root_taxon_specific + '/03_digest/' + goi)
                    shutil.rmtree(root_taxon_specific + '/04_fasta/' + goi)
                    shutil.rmtree(root + 'projects/' + p_id + '/msfragger_secit/' + goi)
                    os.mkdir(root + 'projects/' + p_id + '/msfragger_secit/' + goi + '_did_not_work_retry!')

        else:
            pass

    # print new metadata_species_only to csv
    ncbi_species_df = metadata[metadata[get_accession_assembly_column(metadata)].apply(lambda i: i in species_download.values())]
    len_ncbi_species = len(ncbi_species_df)
    print('f{len_ncbi_species} proteomes in metadata_species file!')
    ncbi_species_df.to_csv(root_taxon_specific + '/00_metafiles/' + goi + '/metadata_genus_specific.csv')
    f_string = len(os.listdir(digest_path))
    print(f'Download finished!')


##################### functions for second iteratino download end ####################################


if __name__ == "__main__":

    # Variables

    print('\n\n#### Action required ####') 	
    root = input('What is your working directory?\nPlease do not use quotation marks.\n')
    

    assembly_path = root + '/00_metafiles/assembly_summary.txt'
    ani_path = root + '/00_metafiles/ANI_report_prokaryotes.txt'

    metadata_output_path = root + '/00_metafiles/'
    proteome_path = root +  '/01_proteomes/'
    unzip_path = root + '/02_unzipped_files/'
    digest_path = root + '/03_digest/'
    fasta_path = root + '/04_fasta/'

    assembly_metadata_path = root + '/00_metafiles/assembly_summary_annot_withANI.csv'
    ## Preprocessing metafile
    metadata_preprocessing(assembly_path, ani_path)

    ## Download files from ncbi
    ncbi_download_reference(assembly_metadata_path , proteome_path, metadata_output_path)

    ## Unzip files
    unzip_files(proteome_path, unzip_path)

    ## digest
    run_digest(0, unzip_path, digest_path)

    ## Create fasta
    write_fasta_2(digest_path, fasta_path)
