

##### MS2Bac by Miriam Abele, refactored by Amirhossein Sakhteman 27 July 2026 #####
1. you need to create the environment The environment is only used for the first step of downloading fasta files 
- Create a conda environment, e.g. conda env create -f linux_conda.yml



2. activate the environment by:
```
conda activate MS2Bac
```

3. Download reference database. This script must be run only once for implementation or when an update is required. 
  - Please respect all terms and conditions from NCBI when downloading data. For further information, see: https://www.ncbi.nlm.nih.gov/home/about/policies/#scripting
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt and store in yourdirectory/00_metafiles
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt and store in yourdirectory/00_metafiles
  - Navigate in Anaconda prompt to your favorite folder\scripts\
  - to check other organisms dynasty https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```
python download_database_MS2Bac.py
```
    Enter your favorite directory when prompted, e.g. C:/Users/MS2Bac_user/Desktop/MS2Bac/. Do not use ' or "


4. Making the mapping database for the peptides and organism
```
python mapping_file_generator.py  # this will make two picke files "idb.pkl.gz" and  "bacci_finder.pkl.gz"   

```
- all the above steps are to make the database and mapping files and needs to be done i.e every year
- in case you need to add a new organism to the current database do as bellow: i.e adding human fasta file

```
python add_new_organism_peptides.py --fasta uniprotkb_proteome_UP000005640_canonical_SwissProt_homo_sapiens_2026_01_09_20417prot.fasta --taxon-id 99999 --genus-tax-id 9606 --organism-name "Homo Sapiens"  --idb idb.pkl.gz  --bacci-finder bacci_finder.pkl.gz --output-fasta peptides_human_fasta.fasta

```

5. by default 4 fasta files are generated in `04_fasta` and normally cause the search to fail, in this case in this folder
```
mkdir search_space
cat *.fasta > /search_space/all_fasta.fasta

```
copy the scripts fragger_run.sh,make_parameter_file.py and chunk_fasta_file.py to search_space folder

```
cd search_space
python chunk_fasta_file all_fasta.fasta 10  # this will make 10 fasta files
rm all_fasta.fasta                          # remove the big fasta from the search space
```

6. Copy the raw file to search_space and run
```
./script.sh <raw_file> 

```
in the search space above the fasta files and raw file should be exisiting, in case you added a new organism accordin to step 4 please add the digested_fasta file also here

7. To combine all the fasta files
```
python combine_fragger_results.py <path_to_psm_files> <name_of_experiment>
```

8. Identification
```
conda deactivate # you need newer version of python to read the zippped picke files
python identification.py <path_to_psm_file from step 7> <path_to_idb_file_from_step 4> <path_to_baccifinder_from_step 4>
```


   