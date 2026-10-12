import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
from statsmodels.stats.multitest import multipletests

"""
usage: python visualize_output.py -i <path_to_output_table_from_identification.py> -o <path_to_outputpdf>

"""


def simulate_relative_drops(
    df,
    count_col="all peptides",
    n_sim=10000,
    seed=42
):
 
    # Keep original dataframe and original row order
    result = df.copy()
 
    # --------------------------------------------------
    # 1. Sort organisms internally by peptide count
    #    but preserve original index
    # --------------------------------------------------
    obs = (
        df.sort_values(count_col, ascending=False)
          .copy()
    )
 
    counts = obs[count_col].to_numpy()
    n_peptides = int(counts.sum())
    n_organisms = len(counts)
 
    # --------------------------------------------------
    # 2. Relative drop compared with FIRST organism
    # --------------------------------------------------
    first_count = counts[0]
 
    obs["relative_drop"] = (
        first_count - obs[count_col]
    ) / first_count
 
    # --------------------------------------------------
    # 3. Null hypothesis
    # Every peptide has equal probability of being
    # assigned to each organism
    # --------------------------------------------------
    p_null = np.repeat(
        1 / n_organisms,
        n_organisms
    )
 
    rng = np.random.default_rng(seed)
 
    # --------------------------------------------------
    # 4. Simulate null distributions
    # --------------------------------------------------
    null_drops = np.zeros(
        (n_sim, n_organisms)
    )
 
    for i in range(n_sim):
 
        sim_counts = rng.multinomial(
            n=n_peptides,
            pvals=p_null
        )
 
        # Rank each simulation
        sim_counts = np.sort(sim_counts)[::-1]
 
        sim_first = sim_counts[0]
 
        # Each rank compared with simulated rank 1
        null_drops[i] = (
            sim_first - sim_counts
        ) / max(sim_first, 1)
 
    # --------------------------------------------------
    # 5. Empirical p-values
    # --------------------------------------------------
    obs["p_value"] = np.nan
 
    # First organism is reference, so start from rank 2
    for rank in range(1, n_organisms):
 
        observed_drop = obs["relative_drop"].iloc[rank]
 
        p = (
            np.sum(
                null_drops[:, rank] >= observed_drop
            ) + 1
        ) / (n_sim + 1)
 
        obs.iloc[
            rank,
            obs.columns.get_loc("p_value")
        ] = p
 
    # --------------------------------------------------
    # 6. Benjamini-Hochberg FDR correction
    # --------------------------------------------------
    obs["FDR"] = np.nan
 
    corrected = multipletests(
        obs["p_value"].iloc[1:],
        method="fdr_bh"
    )[1]
 
    obs.iloc[
        1:,
        obs.columns.get_loc("FDR")
    ] = corrected
 
    # --------------------------------------------------
    # 7. Add results back to ORIGINAL dataframe
    # --------------------------------------------------
    result["relative_drop"] = obs["relative_drop"]
    result["p_value"] = obs["p_value"]
    result["FDR"] = obs["FDR"]
    result['p_FDR'] = -np.log10(result['FDR'].replace(0, np.nan))
    return result, null_drops










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

    df["relative_drop"] = ((df["all peptides"] - df["all peptides"].shift(-1))/ df["all peptides"]).shift(1, fill_value=0)
    truncated_summary, _ = simulate_relative_drops(
    df,
    count_col="all peptides",
    n_sim=10000
)
    fig = plt.figure(figsize=(10,5))
    gs = gridspec.GridSpec(2, 2)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 0])
    
    ax = [ax0, ax1, ax2]
    # getting the plot data in this part

    identified_genus = truncated_summary[~truncated_summary['Organism name'].str.contains('Homo Sapiens')].iloc[0]['Organism name'].split(' ')[0]
    truncated_summary['identified_species'] = truncated_summary['Organism name'].str.contains(identified_genus)
    print(identified_genus)
    sns.scatterplot(data=truncated_summary, x='relative_drop', y="all peptides",ax=ax[0], hue='identified_species')


    for i in range(len(truncated_summary)):
        x = truncated_summary['relative_drop'].iloc[i]
        y = truncated_summary["all peptides"].iloc[i]
        label = truncated_summary["Organism name"].iloc[i]
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
    ax[0].set_xlabel('relative difference to the first organism')



    plot_df = truncated_summary.sort_values("all peptides", ascending=True).tail(20).copy()

    y = np.arange(len(plot_df))

    # Scale p_FDR into visually useful marker sizes
    p = plot_df["p_FDR"].fillna(0).to_numpy()

    if p.max() > p.min():
        circle_size = 5 + 40 * (p - p.min()) / (p.max() - p.min())
    else:
        circle_size = np.full(len(p), 150)




    # Lollipop stems
    ax[1].hlines(
        y=y,
        xmin=0,
        xmax=plot_df["all peptides"],
        color="lightgray",
        linewidth=2
    )

    # Lollipop circles
    sc = ax[1].scatter(
        plot_df["all peptides"],
        y,
        s=circle_size,
        c=plot_df["p_FDR"],
        cmap="viridis",
        edgecolor="black",
        linewidth=0.6,
        zorder=3
    )

    # Organism labels
    ax[1].set_yticks(y)
    ax[1].set_yticklabels(
        plot_df["Organism name"],
        fontsize=5
    )

    # Keep y-axis tightly around ALL data
    ax[1].set_ylim(-0.5, len(plot_df) - 0.5)

    ax[1].set_xlabel("All peptides")
    ax[1].set_ylabel("")
    ax[1].set_title("Peptide assignments by organism")

    # Colour bar
    cbar = plt.colorbar(
        sc,
        ax=ax[1],
        pad=0.02
    )
    cbar.set_label("-log10(FDR)")

    # Remove unnecessary borders
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["right"].set_visible(False)

    #print(truncated_summary)
    plot_df = truncated_summary.sort_values(by='all peptides', ascending=False)[:10][['all peptides',  'Organism name','FDR']].copy()
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

