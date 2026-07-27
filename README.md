

##### MS2Bac by Miriam Abele #####
1. you need to create the environment The environment is only used for the first step of downloading fasta files 
```
  - Create a conda environment, e.g. conda env create -f linux_conda.yml
```


2. activate the environment by:
```
conda activate MS2Bac
```

3. Download reference database. This script must be run only once for implementation or when an update is required. 
  - Please respect all terms and conditions from NCBI when downloading data. For further information, see: https://www.ncbi.nlm.nih.gov/home/about/policies/#scripting
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt and store in yourdirectory/00_metafiles
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt and store in yourdirectory/00_metafiles
  - Navigate in Anaconda prompt to your favorite folder\scripts\
  - Start script: python download_database_MS2Bac.py

    Enter your favorite directory when prompted, e.g. C:/Users/MS2Bac_user/Desktop/MS2Bac/. Do not use ' or "


4. Making the mapping database for the peptides and organism
```
python mapping_file_generator.py  # this will make two picke files

```


5. by default 4 fasta files are generated and normally cause the search to fail, in this case
```
mkdir search_space
cat *.fasta > /search_space/all_fasta.fasta

```
copy the scripts fragger_run.sh,make_parameter_file.py and chunk_fasta_file.py to search_space folder

```
python chunk_fasta_file all_fasta.fasta 10  # this will make 10 fasta files
rm all_fasta.fasta                          # remove the big fasta from the search space
```
6. Copy the raw file to search_space and run
```
./script.sh BBM_428_P110_31_MIA_007_30.raw # replce your fasta file, this will run the searches using Fragpipe 24 and makes the PSM files per chunk

```
7. To combine all the fasta files
```
python combine_fragger_results.py <path_to_psm_files> <name_of_experiment>
```

8. Identification
```
conda deactivate # you need newer version of python to reada the zippped picke files

```


   