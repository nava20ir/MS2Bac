# pipeline to run MS2BAC with casanovo

## it does the search using casanovo
## Instead of search against database of digested proteins, it counts the number of times a pepeide matches to any of the proteins of that organism

you need one file named `all_proteome.csv` with three columns 
`Proteins` is the FASTA header
`sequence`  full sequence of the protein  
`Organism` the organism name refering to that protein

- To run the pipeline

```
./batch_convert.sh # this will convert all raw files into mzml
# use your casanovo environment
conda activate casanovo
./batch_casanovo.sh
./batch_prepare.sh  # this will do the final report creation
```

