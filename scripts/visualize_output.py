import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np

"""
usage: python visualize_output.py -i <path_to_output_table_from_identification.py> -o <path_to_outputpdf>

"""

def main(output_summary, output_pdf):
    '''
    Function visualizes bacterium identification in a swarmplot (one for unique peptides and one for all)

    Parameters
    ----------
    output_summary : pandas DataFrame
        file with all bacteria and raw files and the corresponding peptide number and other data

    output_pdf: pdf object from packages PdfPages

    Returns
    --------
        Writes pdf files
    '''

    identification = pd.read_csv(output_summary)
    #identification = identification[~identification['Organism name'].str.contains('Homo Sapiens')]
    with PdfPages(output_pdf) as pdf:
        visualize_identification(identification.reset_index().sort_values(by='all peptides', ascending=False), pdf)





def visualize_identification(summary, pdf):
    '''
    Function visualizes bacterium identification in a swarmplot (one for unique peptides and one for all)

    Parameters
    ----------
    summary : pandas DataFrame
        file with all bacteria and raw files and the corresponding peptide number and other data

    pdf: pdf object from packages PdfPages

    Returns
    --------
        Writes pdf files
    '''

    for e, g in summary.groupby('Experiment'):

        fig = plt.figure(figsize=(10,5))
        gs = gridspec.GridSpec(2, 2)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[1, :])
        ax = [ax0, ax1, ax2]

        
        truncated_summary = summary[summary['Experiment'] == e].copy()
        identified_genus = truncated_summary[~truncated_summary['Organism name'].str.contains('Homo Sapiens')].iloc[0]['Genus']
        truncated_summary['identified_species'] = (truncated_summary['Genus'] == identified_genus)

        sns.stripplot(data=truncated_summary[~truncated_summary['Organism name'].str.contains('Homo Sapiens')], x='Experiment', y='all peptides', ax=ax[0], hue='identified_species')
        ax[0].set_title('all peptides')
        ax[0].set_xlabel('')

        sns.stripplot(data=truncated_summary[~truncated_summary['Organism name'].str.contains('Homo Sapiens')], x='Experiment',y='unique peptides', ax=ax[1], hue='identified_species')
        ax[1].set_title('unique sequences')
        ax[1].set_xlabel('')
        #print(truncated_summary)
        plot_df = truncated_summary.sort_values(by='all peptides', ascending=False)[:10][['score', 'all peptides', 'relative', 'unique peptides', 'Genus', 'Organism name']].copy()
        filtered_df = plot_df[~plot_df['Organism name'].str.contains('Homo Sapiens')]
        print(plot_df[plot_df['Organism name'].str.contains('Homo Sapiens')])
        ccol = plt.cm.BuPu(np.full(len(filtered_df.columns), 0.1))
        table = ax[2].table(cellText=plot_df.values,
                        colLabels=plot_df.columns,
                        loc='center',
                        colColours=ccol)
        table.auto_set_column_width(col=list(range(len(plot_df.columns))))
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1,1)

        ax[2].set_frame_on(False)
        ax[2].axes.get_yaxis().set_visible(False)
        ax[2].axes.get_xaxis().set_visible(False)

        plt.tight_layout()
        
        pdf.savefig(fig)
        plt.close()
    print('Visualization is done. ')


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Visualize bacterium identification in a swarmplot (one for unique peptides and one for all)')
    parser.add_argument('-i', '--input_summary', type=str, required=True, help='Input summary file with all bacteria and raw files and the corresponding peptide number and other data')
    parser.add_argument('-o', '--output_pdf', type=str, required=True, help='Output pdf file to save the plots')
    args = parser.parse_args()

    main(args.input_summary, args.output_pdf)
