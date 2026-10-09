import ahocorasick
import pandas as pd

from collections import Counter

def main(peptide_file):
    peptides_df = pd.read_csv(peptide_file)
    print('reading frames')
    proteins_df = pd.read_csv('all_proteome.csv')
    print('Runnig code')
    proteins_df['sequence'] = proteins_df['sequence'].str.replace('I|L','J',regex=True)
    # Get unique peptide sequences
    peptides = peptides_df["Sequence"].dropna().unique()
    # Create Aho-Corasick automaton
    A = ahocorasick.Automaton()

    for i, peptide in enumerate(peptides):
        A.add_word(peptide, i)

    A.make_automaton()

    # Count peptide matches for each organism
    counts = Counter()

    for protein, organism in zip(
        proteins_df["sequence"].values,
        proteins_df["Organism"].values
    ):
        matches = set()

        # Find all unique peptides present in this protein
        for _, peptide_id in A.iter(protein):
            matches.add(peptide_id)

        # Each peptide is counted only once per protein
        counts[organism] += len(matches)

    # Print results
    count_df = pd.DataFrame.from_dict(counts, orient='index', columns=['Peptide_Count'])
    count_df = count_df.sort_values(by='Peptide_Count', ascending=False)
    #count_df["relative_drop"] = (count_df["Peptide_Count"] - count_df["Peptide_Count"].shift(-1)) / count_df['Peptide_Count']
    count_df.to_csv(f'count_per_organism_protein_level_{peptide_file}')


if __name__ == "__main__":
    import sys
    peptide_file = sys.argv[1]
    main(peptide_file)
