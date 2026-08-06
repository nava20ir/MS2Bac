# MS2Bac

MS2Bac identifies bacterial and fungal species from MS/MS peptide spectra using an MSFragger-based search workflow.

Original implementation by Miriam Abele and Christina Ludwig.

Refactored and extended by Amirhossein Sakhteman (July 2026).

---

# Installation

## 1. Create the Conda Environment

This environment is required for database generation and maintenance.

```bash
conda env create -f linux_conda.yml
```

Activate the environment:

```bash
conda activate MS2Bac
```

---

## 2. Pull the MSFragger Docker Image

MSFragger is executed via Docker.

```bash
docker pull nava20ir/fragpipebbm:latest
```

---

# Directory Structure

```text
MS2Bac/
├── Makefile
├── linux_conda.yml
├── 00_metafiles/
├── 01_proteomes/
├── 02_unzip/
├── 03_digest/
├── 04_fasta/
├── results/
│
├── download_database_MS2Bac.py
├── mapping_file_generator.py
├── add_new_organism_peptides.py
├── chunk_fasta_file.py
├── make_fragger_parameter_file.py
├── combine_fragger_results.py
├── identification.py
└── visualize_output.py
```

---

# Quick Start

Once the database has been generated, process a RAW file using:

```bash
make report \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample \
    IDB=/path/to/idb.pkl.gz \
    BACCI=/path/to/bacci_finder.pkl.gz
```

The final report will be generated in:

```text
results/sample/sample_identification.pdf
```

---

# Database Generation

This step only needs to be run:

- during initial installation
- when updating the reference database
- when rebuilding the database for another organism group

---

## 1. Download Metadata Files

Create:

```text
00_metafiles/
```

Download the following files from NCBI.

### Bacteria

```text
assembly_summary.txt
ANI_report_prokaryotes.txt
```

### Fungi (Saccharomycotina)

```text
assembly_summary.txt
ANI_report_saccharomycotina.txt
```

Reference:

https://ftp.ncbi.nlm.nih.gov/genomes/refseq/

Please respect all NCBI download policies.

---

## 2. Download and Build the Reference Database

### Bacteria

```bash
make database \
    ANI_FILE=00_metafiles/ANI_report_prokaryotes.txt
```

### Fungi

```bash
make database \
    ANI_FILE=00_metafiles/ANI_report_saccharomycotina.txt
```

This step:

- downloads representative proteomes
- unzips FASTA files
- performs in-silico digestion
- creates peptide FASTA databases

---

## 3. Generate Mapping Databases

```bash
make mapping
```

This generates:

```text
idb.pkl.gz
bacci_finder.pkl.gz
```

These files are required for organism identification.

---

# Adding a New Organism

Additional proteomes can be added to an existing database.

Example:

```bash
python add_new_organism_peptides.py \
    --fasta human.fasta \
    --taxon-id 99999 \
    --genus-tax-id 9606 \
    --organism-name "Homo sapiens" \
    --idb idb.pkl.gz \
    --bacci-finder bacci_finder.pkl.gz \
    --output-fasta peptides_human.fasta
```

This will:

1. Digest the FASTA file.
2. Update `idb.pkl.gz`.
3. Update `bacci_finder.pkl.gz`.
4. Generate an optional peptide FASTA file.

---

# Preparing the Search Space

Create FASTA chunks for MSFragger searching:

```bash
make prepare-search
```

or customize the number of chunks:

```bash
make prepare-search N_CHUNKS=20
```

This step:

1. Combines FASTA files from `04_fasta/`
2. Creates `04_fasta/search_space/`
3. Splits the search space into chunks
4. Prepares FASTA files for MSFragger

---

# Running an Identification

Run the complete workflow:

```bash
make report \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample \
    IDB=idb.pkl.gz \
    BACCI=bacci_finder.pkl.gz
```

Pipeline:

```text
RAW
 ↓
MSFragger Search
 ↓
combine_fragger_results.py
 ↓
psm.csv
 ↓
identification.py
 ↓
output_table.csv
 ↓
visualize_output.py
 ↓
PDF Report
```

---

# Running Individual Pipeline Steps

## Search Only

```bash
make search \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample
```

---

## Combine Search Results

```bash
make combine \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample
```

---

## Identification Only

```bash
make identify \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample \
    IDB=idb.pkl.gz \
    BACCI=bacci_finder.pkl.gz
```

---

## Generate PDF Only

```bash
make visualise EXPERIMENT=sample
```

---

# Useful Commands

Display all available Make targets:

```bash
make help
```

Remove generated FASTA search space:

```bash
make clean-search
```

Remove results from one experiment:

```bash
make clean-results EXPERIMENT=sample
```

---

# Output Files

For experiment `sample`:

```text
results/sample/
├── psm.csv
├── output_table.csv
├── summary_output.csv
└── sample_identification.pdf
```

---

# Advanced Options

Change the number of FASTA chunks:

```bash
make report \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample \
    N_CHUNKS=20 \
    IDB=idb.pkl.gz \
    BACCI=bacci_finder.pkl.gz
```

Adjust Java memory for MSFragger:

```bash
make report \
    RAW=/path/to/sample.raw \
    EXPERIMENT=sample \
    JAVA_MEM=64G \
    IDB=idb.pkl.gz \
    BACCI=bacci_finder.pkl.gz
```

---

# Notes

- Database generation typically only needs to be performed once.
- Identification can be run repeatedly using the generated `idb.pkl.gz` and `bacci_finder.pkl.gz`.
- MSFragger is executed through Docker (`nava20ir/fragpipebbm:latest`).
- Results from each 

##### MS2Bac by Miriam Abele, refactored by Amirhossein Sakhteman 27 July 2026 #####
1. you need to create the environment The environment is only used for the first step of downloading fasta files 
- Create a conda environment, e.g. conda env create -f linux_conda.yml



2. activate the environment by:
```
conda activate MS2Bac
```
- for the msfragger search we use the dockerized image in the docker hub
- first pull the image; this image contains all files needed to run fragpipe
```
docker pull nava20ir/fragpipebbm:latest
```

3. Download reference database. This script must be run only once for implementation or when an update is required. 
  - Please respect all terms and conditions from NCBI when downloading data. For further information, see: https://www.ncbi.nlm.nih.gov/home/about/policies/#scripting
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt and store in yourdirectory/00_metafiles
  - Download: https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt and store in yourdirectory/00_metafiles
  - Navigate in Anaconda prompt to your favorite folder\scripts\
  - to check other organisms dynasty https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
```
python download_database_MS2Bac.py <path-to-folder> <path-to-ANI file for that organism>
python download_database_MS2Bac.py /home/asakhteman/ms2bac/ ANI_report_saccharomycotina.txt # this is for fungus
python download_database_MS2Bac.py /home/asakhteman/ms2bac/ ANI_report_prokaryotes.txt # this is for Bacteria
```
   


4. Making the mapping database for the peptides and organism
```
python mapping_file_generator.py <path-to-folder>   # this will make two picke files "idb.pkl.gz" and  "bacci_finder.pkl.gz"   

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


6. Copy the raw file to search_space and run this command once per raw file
```
./script.sh <raw_file> 

```
in the search space above the fasta files and raw file should be exisiting, in case you added a new organism accordin to step 4 please add the digested_fasta file also here

7. To combine all the generated PSM  files in to one psm.file
```
python combine_fragger_results.py <path_to_psm_files> <name_of_experiment> # this will make the final psm file
```

- Using Casanovo;
- alternatively one can use casanovo to get the PSMs in this case for the step 6 and 7 use the bellow codes

```
conda activate casanovo
python raw2mzml.py -i [PATH_TO MZML] -o [output]
cd output 
casanovo sequence [PATH_TO MZML]

```



8. Identification
```
conda deactivate                             # you need newer version of python to read the zippped picke files not the old env
python identification.py <path_to_psm_file from step 7> <path_to_idb_file_from_step 4> <path_to_baccifinder_from_step 4>
```

9. visualization
```
python visualize_output.py -i <path_to_output_table_from_identification.py> -o <path_to_outputpdf>

```


   