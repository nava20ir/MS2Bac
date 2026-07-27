

##### MS2Bac by Miriam Abele #####
  1. you need to create the environment
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
	### Please note ### 
	Follow NCBI instructions for download times. See: https://www.ncbi.nlm.nih.gov/home/about/policies/#scripting

    - Crate four fasta.pepindex files by following the following instructions:
       Open fragpipe, add one raw file of choice, change database to first fasta file, set an output directory in the 'Run' tab, then start a "search". Use default parameters. This procedure generates 
       a .pepindex file which is needed to run MSfragger with the command line. Additionally, it will create a params file. 
       Note: MSFragger will notify you that no decoys are included. This is fine. Repeat the procedure for the other three fasta files.

    - Optional: If you want to use another MSfragger version than the tested ones, please move the params file from the MSFragger output directory to C:/Users/MS2Bac_user/Desktop/MS2Bac/05_params/. You must create four params file copies. Please name them fragger0.params, fragger1.params, fragger2.params, fragger3.params
      Please note: You MUST change the parameter "output_format" to "tsv_pepXML", otherwise the script will not work. All other parameters can be adjusted according
      to need. An example can be found here: C:/Users/MS2Bac_user/Desktop/MS2Bac/05_params

     - Now you MUST change the database name in the params file which are located in C:/Users/MS2Bac_user/Desktop/MS2Bac/05_params (either you have created them yourself or you use the default ones). Please edit: database_name = yourdirectory/yyyy-mm-dd_x.fasta with x = 0, 1, 2, 3 (same names as in C:/Users/MS2Bac_user/Desktop/MS2Bac/04_fasta/).
 
    - If you want to use the same version as tested, just make sure that you edit the params file provided in the 05_params folder downloaded from ZENODO. 
      In each of the four params files, you MUST change the database_name = yourdirectory/yyyy-mm-dd_x.fasta with x = 0, 1, 2, 3 (same names as in C:/Users/MS2Bac_user/Desktop/MS2Bac/04_fasta/

4. by default 4 fasta files are generated and normally cause the search to fail, in this case
```
mkdir search_space
cat *.fasta > /search_space/all_fasta.fasta

```
copy the scripts fragger_run.sh,make_parameter_file.py and chunk_fasta_file.py to search_space folder

```
python chunk_fasta_file all_fasta.fasta 10  # this will make 10 fasta files
rm all_fasta.fasta                          # remove the big fasta from the search space
```
5. Copy the raw file to search_space and run
```
./script.sh BBM_428_P110_31_MIA_007_30.raw # replce your fasta file, this will run the searches using Fragpipe 24 and makes the PSM files per chunk

```
6. To combine all the fasta files
```
python combine_fragger_results.py <path_to_psm_files> <name_of_experiment>
```
7. Making the mapping database for the peptides and organism
```
python mapping_file_generator.py  # this will make two picke files

```
8. Identification
```


```


   