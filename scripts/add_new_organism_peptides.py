#!/usr/bin/env python3
"""
Usage:
python add_new_organism_peptides.py --fasta new_organism.faa --rule trypsin --taxon-id 123456 --genus-tax-id 9876 --organism-name "My bacterium"  --idb idb.pkl.gz  --bacci-finder bacci_finder.pkl.gz
"""

import argparse
import gzip
import pickle
import re
from pyteomics import fasta
import pyteomics
from pyteomics.parser import cleave, expasy_rules
from typing import Any, Dict, Iterable, Iterator, List, Optional, Set, Tuple





def digest_fasta(
    fasta_file,
    output_fasta=None,
    rule="trypsin",
    missed_cleavages=2,
    min_length=7,
    max_length=30,
):
    if rule == 'D_cut': digest_rule='((?=D)|(?<=D))'
    else:
        digest_rule = expasy_rules[rule]
    
    peptides = set()

    with open(fasta_file, "rt") as f:
        for protein in pyteomics.fasta.read(f):

            for peptide in cleave(
                protein.sequence,
                digest_rule,
                missed_cleavages,
                min_length,
            ):
                if rule == 'D_cut':
                    peptide = peptide.strip('D')

                if (
                    len(peptide) <= max_length
                    and "X" not in peptide
                    and "U" not in peptide
                    and "B" not in peptide
                    and "Z" not in peptide
                    and "i" not in peptide
                ):
                    peptides.add(peptide)

    if output_fasta is not None:
        fasta_entries = [
            (f"peptide_{i}", peptide)
            for i, peptide in enumerate(sorted(peptides))
        ]
        fasta.write(fasta_entries, output_fasta)

    return peptides


def add_organism_to_idb(idb, organism_peptides, taxon_id):

    for seq in organism_peptides:
        seq = seq.replace("I", "J").replace("L", "J")
        #idb[seq].add(taxon_id)
        idb.setdefault(seq, set()).add(taxon_id)

    return idb


def add_organism_to_bacci_finder(
    bacci_finder,
    taxon_id,
    organism_name,
    genus_tax_id,
    organism_peptides,
):

    bacci_finder[taxon_id] = (
        organism_name,
        genus_tax_id,
        len(set(organism_peptides)),
    )

    return bacci_finder


def main():

    parser = argparse.ArgumentParser(
        description="Add a new organism FASTA to an existing BacciFinder database."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "--taxon-id",
        required=True,
        type=int,
        help="Species taxon ID"
    )

    parser.add_argument(
        "--genus-tax-id",
        required=True,
        type=int,
        help="Genus taxon ID"
    )

    parser.add_argument(
        "--organism-name",
        required=True,
        help="Organism name"
    )

    parser.add_argument(
        "--idb",
        required=True,
        help="idb.pkl.gz"
    )

    parser.add_argument(
        "--rule",
        default="trypsin",
        help="Digestion rule (trypsin, chymotrypsin, D_cut, etc.)"
    )

    parser.add_argument(
        "--bacci-finder",
        required=True,
        help="bacci_finder.pkl.gz"
    )

    parser.add_argument(
        "--output-fasta",
        default=None,
        help="Optional peptide FASTA output"
    )

    args = parser.parse_args()

    print("Digesting FASTA...")

    peptides = digest_fasta(
        args.fasta,
        rule=args.rule,
        output_fasta=args.output_fasta,
    )

    print("Loading bacci_finder...")
    try:
        with gzip.open(args.bacci_finder, "rb") as f:
            bacci_finder = pickle.load(f)
    except FileNotFoundError:
        bacci_finder = {}
        print("bacci_finder file not found. Creating a new one.")

    if args.taxon_id in bacci_finder:
        organism_name = bacci_finder[args.taxon_id][0]

        raise SystemExit(
            f"ERROR: taxon_id {args.taxon_id} already exists "
            f"({organism_name}). Aborting."
        )

    print("Loading idb...")
    try:
        with gzip.open(args.idb, "rb") as f:
            idb = pickle.load(f)
    except FileNotFoundError:
        idb = {}
        print("idb file not found. Creating a new one.")



    print(f"Generated {len(peptides):,} unique peptides")

    print("Updating idb...")

    add_organism_to_idb(
        idb=idb,
        organism_peptides=peptides,
        taxon_id=args.taxon_id,
    )

    print("Updating bacci_finder...")

    add_organism_to_bacci_finder(
        bacci_finder=bacci_finder,
        taxon_id=args.taxon_id,
        organism_name=args.organism_name,
        genus_tax_id=args.genus_tax_id,
        organism_peptides=peptides,
    )

    print("Saving idb...")

    with gzip.open(args.idb, "wb", compresslevel=1) as f:
        pickle.dump(
            idb,
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    print("Saving bacci_finder...")

    with gzip.open(args.bacci_finder, "wb", compresslevel=1) as f:
        pickle.dump(
            bacci_finder,
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    print()
    print("Finished")
    print(f"Taxon ID: {args.taxon_id}")
    print(f"Organism: {args.organism_name}")
    print(f"Peptides: {len(peptides):,}")
    print(f"Total database peptides: {len(idb):,}")
    print(f"Total organisms: {len(bacci_finder):,}")


if __name__ == "__main__":
    main()