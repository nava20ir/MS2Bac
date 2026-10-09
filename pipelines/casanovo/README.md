# pipeline to run MS2BAC with casanovo

you need one file named `all_proteome.csv` with three columns 
`Proteins` is the FASTA header
`sequence`  full sequence of the protein  
`Organism` the organism name refering to that protein

- To run the pipeline

```
./batch_convert.sh # this will convert all raw files into mzml
conda activate casanovo
./batch_casanovo.sh
./batch_prepare.sh  # this will do the final report creation
```

