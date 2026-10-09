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
    
    identification = identification.rename(columns={"Unnamed: 0": "Organism name","Peptide_Count": "all peptides"})
    print(identification.columns)



    #identification = identification[~identification['Organism name'].str.contains('Homo Sapiens')]
    with PdfPages(output_pdf) as pdf:
        visualize_identification(identification.reset_index().sort_values(by='all peptides', ascending=False), pdf)


def visualize_identification(df, pdf):
    '''
    Function visualizes bacterium identification in a swarmplot (one for unique peptides and one for all)

    Parameters
    ----------
    df : pandas DataFrame
        file with all bacteria and raw files and the corresponding peptide number and other data

    pdf: pdf object from packages PdfPages

    Returns
    --------
        Writes pdf files
   '''
    df['Experiment'] = 'test'
    print(df.columns)
    df["relative_drop"] = ((df["all peptides"] - df["all peptides"].shift(-1))/ df["all peptides"]).shift(1, fill_value=0)

    for e, g in df.groupby('Experiment'):

        fig = plt.figure(figsize=(10,5))
        gs = gridspec.GridSpec(2, 1)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[1, 0])
        
        ax = [ax0, ax1]
        # getting the plot data in this part
        truncated_summary = df[df['Experiment'] == e].copy()
        identified_genus = truncated_summary[~truncated_summary['Organism name'].str.contains('Homo Sapiens')].iloc[0]['Organism name'].split(' ')[0]
        truncated_summary['identified_species'] = truncated_summary['Organism name'].str.contains(identified_genus)
        print(identified_genus)
        sns.scatterplot(data=truncated_summary, x='relative_drop', y="all peptides",ax=ax[0], hue='identified_species')


        for i in range(len(df)):
            x = df['relative_drop'].iloc[i]
            y = df["all peptides"].iloc[i]
            label = df["Organism name"].iloc[i]
            if i < 1:  # Annotate only the first 10 points
                ax[0].annotate(
                    label,
                    xy=(x, y),                 # point to annotate
                    xytext=(x + 0.0001 , y + 1 ),   # label position
                    arrowprops=dict(
                        arrowstyle="->",
                        lw=0.5
                    ),
                    fontsize=6
                )
        ax[0].set_title('all peptides')
        ax[0].set_xlabel('')
        #print(truncated_summary)
        plot_df = truncated_summary.sort_values(by='all peptides', ascending=False)[:10][['all peptides',  'Organism name']].copy()
        filtered_df = plot_df[~plot_df['Organism name'].str.contains('Homo Sapiens')]
        print(plot_df[plot_df['Organism name'].str.contains('Homo Sapiens')])
        ccol = plt.cm.BuPu(np.full(len(filtered_df.columns), 0.1))
        table = ax[1].table(cellText=plot_df.values,
                        colLabels=plot_df.columns,
                        loc='center',
                        colColours=ccol)
        table.auto_set_column_width(col=list(range(len(plot_df.columns))))
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1,1)

        ax[1].set_frame_on(False)
        ax[1].axes.get_yaxis().set_visible(False)
        ax[1].axes.get_xaxis().set_visible(False)

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
