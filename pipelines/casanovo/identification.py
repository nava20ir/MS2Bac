import pickle
import gzip
from collections import defaultdict
import pandas as pd




def count_seqs(experimental_seqs, database, bacci_finder):

    # Assign variables
    bacci_score_counter = defaultdict(lambda: 0)
    bacci_unique_counter = defaultdict(lambda: 0)
    bacci_total_counter = defaultdict(lambda: 0)
    bacci_psms_counter = defaultdict(lambda: 0)


    # iterate through all identified sequences from experiment and count
    for seq, score  in experimental_seqs:
        print(seq, score)
        try:
            for organism_id in database[seq]:
                
                bacci_score_counter[organism_id] += score
                bacci_total_counter[organism_id] += 1
                bacci_psms_counter[organism_id] = len(experimental_seqs)

            # counts only unique peptides (unique if database value equals 1)
            if (len(database[seq]) == 1):
                bacci_unique_counter[organism_id] += 1
        except KeyError:
            print(f"Sequence {seq} not found in database. Skipping.")    

    # Add missing organism_ids and give them a score = 0, a count = 0, the correct number of total psms, and unique = 0 
    for i in bacci_finder.keys():
        bacci_score_counter.setdefault(i, 0)
        bacci_total_counter.setdefault(i, 0)
        bacci_psms_counter.setdefault(i, len(experimental_seqs))
        bacci_unique_counter.setdefault(i, 0)

    # convert defaultdict into pandas Dataframe and sort them ascending
    sc = pd.DataFrame(pd.Series(bacci_score_counter), columns=['score'])
    uc = pd.DataFrame(pd.Series(bacci_unique_counter), columns=['unique peptides'])
    tc = pd.DataFrame(pd.Series(bacci_total_counter), columns=['all peptides'])
    psm = pd.DataFrame(pd.Series(bacci_psms_counter), columns=['psms'])
    

    # merge all counts and unique counts to one dataframe and return
    sc = pd.merge(left=sc, right=uc, left_index=True, how='outer', right_index=True)
    sc = pd.merge(left=sc, right=psm, left_index=True, how='outer', right_index=True)
    final_df = pd.merge(left=sc, right=tc, left_index=True, how='outer', right_index=True)
    print(final_df)
    return final_df


def summary_output(summary, bacci_finder):
    try:

        index = pd.Index(((e, df_index) for e, df in summary.items() for df_index in df.index),
                        name=('Experiment', 'Organism'))

        score = pd.concat([df['score'] for e, df in summary.items()],
                        ignore_index=True)

        unique_counts = pd.concat([df['unique peptides'] for e, df in summary.items()],
                                ignore_index=True)

        total = pd.concat([df['all peptides'] for e, df in summary.items()],
                                ignore_index=True)

        psms = pd.concat([df['psms'] for e, df in summary.items()],
                                ignore_index=True)


        df = pd.DataFrame({'score':score.values, 'all peptides': total.values, 'psms': psms.values, 'unique peptides':unique_counts.values}, index=index).reset_index()

        df['map'] = df['Organism'].map(bacci_finder)
        df[['Organism name','Genus', '# peptides']] = pd.DataFrame(df['map'].tolist(),index=df.index)
        df['relative'] = df['all peptides'] / df['# peptides']
        df.to_csv('output_table.csv')
        #path_test = os.path.join(output_path, 'output_table.csv')

        #print(f'Outputtable printed to {path_test}')
        short = df.groupby('Experiment').apply(lambda g: g.sort_values(by='all peptides', ascending=False).iloc[0])
        return short
    except ValueError:
        print('No sequences were identified. Please check your input file.')
        return None








if __name__ == "__main__":

    import sys
    psm_file = sys.argv[1]
    idb_file = sys.argv[2]
    bacci_finder_file = sys.argv[3]
    psms = pd.read_csv(psm_file)
    print('reading bacci finder')
    with gzip.open(bacci_finder_file, "rb") as f: bacci_finder = pickle.load(f)
    print('reading big file')
    with gzip.open(idb_file, "rb") as f: idb = pickle.load(f)
    print('going with final step')
    # Prepare psms file
    d = defaultdict(list)
    psms = psms.groupby(['Experiment', 'Sequence'], as_index=False)['Score'].max()
    for tup in psms.itertuples():  d[tup.Experiment].append((tup.Sequence, tup.Score))
    sample_seqs = pd.Series(dict(d))
    res = summary_output(sample_seqs.apply(count_seqs, database=idb, bacci_finder=bacci_finder), bacci_finder=bacci_finder)
    res.to_csv('summary_output.csv', index=False)
