def add_organism_to_idb(idb, organism_peptides, taxon_id):
    """
    Add a new organism's peptides to an existing idb.

    Parameters
    ----------
    idb : defaultdict(set)
        Peptide -> set(taxon_ids) database.
    organism_peptides : iterable
        Collection of peptide sequences.
    taxon_id : int or str
        Taxon ID of the organism.

    USAGE:
    new_taxon_id = 12345
    organism_peptides = ["PEPTIDE", "ABCDEF", "XYZ123"]
    idb = add_organism_to_idb(idb, organism_peptides, new_taxon_id)
    """
    for seq in organism_peptides:
        seq = seq.replace("I", "J").replace("L", "J")
        idb[seq].add(taxon_id)

    return idb



if __name__ == '__main__':
    print('reading bacci_finder filer')
    with gzip.open(bacci_finder_file, "rb") as f: bacci_finder = pickle.load(f)

    print('reading database picle file')
    with gzip.open(idb_file, "rb") as f: idb = pickle.load(f)


