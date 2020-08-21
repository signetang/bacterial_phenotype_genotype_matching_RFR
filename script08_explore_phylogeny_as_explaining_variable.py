import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
import pandas as pd
import seaborn as sns
from script02_color_scheme import load_color_scheme

sns.set()


def make_dendrogram_and_clustermap(genomic_matrix, colors, genomic_type, output, target, phenotypes):

    # Making a vector of CHCC numbers and their corresponding colors (colored by sub species)
    taxonomy = pd.DataFrame() # save the taxonomy (indexed by CHCC numbers) separate from the dataframe 'target'
    taxonomy["Taxonomy"] = genomic_matrix["Taxonomy"]
    row_colors = taxonomy['Taxonomy'].map(colors)

    # Plot the dendrogram of the genomic_matrix and the clustermap of the genomic_matrix, and getting the "phylogenetic" row order:
    print 'Making dendrogram ("phylogenetic tree")...'

    fig = plt.figure(figsize=(16,int(len(genomic_matrix)/10.)))
    dist = hierarchy.distance.pdist(genomic_matrix.iloc[:,1:]) #exclude the taxonomy column of the genomic_matrix
    link = hierarchy.linkage(dist, method='weighted', metric='euclidean')
    dend = hierarchy.dendrogram(link, labels=genomic_matrix.index, color_threshold=15, above_threshold_color='grey', orientation='left', leaf_font_size=15)
    ax = plt.gca()
    # Display the strain name in the relevant color:
    strain_labels = ax.get_ymajorticklabels()
    for c in strain_labels:
        c.set_color(row_colors[c.get_text()])
        c.set_fontsize(9)
    ax.set_title('Phylogenetic tree based on '+genomic_type, fontsize=35)
    ax.invert_yaxis() # for some reason the dendrogram is oposite the other plots
    plt.tight_layout()
    fig.savefig(output+'_dendrogram.png')

    print 'Making clustermap...'
    plt.figure()
    cm = sns.clustermap(genomic_matrix.iloc[:,1:], row_linkage=link, figsize=(int(len(genomic_matrix)/3.),int(len(genomic_matrix)/8.)), yticklabels=True, cmap="Greys") #, row_colors=row_colors)
    strain_labels = cm.ax_heatmap.get_yticklabels()
    for c in strain_labels:
        c.set_color(row_colors[c.get_text()])
        c.set_fontsize(9)
    cm.savefig(output+'_clustermap.png')


    # Okay, try to keep up here.. Lots of reorderings..
    # The issue is that the rows of target are not in the same order as the rows of genomic_matrix.
    # Therefore the reindexing to the clustermap order (the "phylogenetic" row order) will land them completely differently!

    # Getting the target row order to be the same as the genomic_matrix row order, so that it can be reordered by the "phylogenetic" reordering:
    target = target.reindex(genomic_matrix.index) # get the order of target to where the genomic_matrix was before the clustermap clustering
    target['Strains'] = target.index # make copy of Strains column before overwriting it with an integer index
    target.index = np.arange(len(target)) # make integer index
    target = target.reindex(cm.dendrogram_row.reordered_ind) # apply the same reordering that clustermap did to genomic_matrix
    target = target.set_index('Strains') # make the Strains column the index again
    target['Taxonomy'] = genomic_matrix['Taxonomy']

    row_colors = row_colors.rename('') # Avoid plotting 'Taxonomy' under the row colors

    # Plotting the clustermap of the phenotypes (without clustering, as the row order is already correct)
    for pheno in phenotypes:
        print 'Making '+pheno+' clustermap now...'
        columns = [column for column in target.columns if column[:len(pheno)]==pheno]

        plt.figure()
        cm = sns.clustermap(target[columns], row_colors=row_colors, row_cluster=False, col_cluster=False, figsize=(40,90))
        cm.cax.set_position((0.15, 0.8345, 0.025, 0.0425)) # colorbar width adjusted
        plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), rotation=0) # x-labels horizontal
        cm.ax_heatmap.xaxis.set_label_position('top') # Adds the space for the labels on top
        cm.ax_heatmap.yaxis.tick_left() # Puts the labels on the left side
        cm.ax_heatmap.xaxis.tick_top() # Puts the labels on top
        col_ax = cm.ax_col_dendrogram
        sns.boxplot(data=target[columns], ax=col_ax, color='azure')

        # Put the "taxonomic" dendrogram in the plot:
        row_ax = cm.ax_row_dendrogram
        dend = hierarchy.dendrogram(link, color_threshold=9, above_threshold_color='grey', orientation='left', ax=row_ax)
        row_ax.invert_yaxis() # for some reason the dendrogram is oposite the other plots

        # Make all the phenotype values go in one column and a description of the phenotype in another column:
        target['Strains'] = target.index
        target2 = pd.melt(target, id_vars=['Strains', 'Taxonomy'], var_name="Variables", value_name="Value")

        columns2 = [column for column in target.columns if column[:len(pheno)]!=pheno]
        target2 = target2[~target2['Variables'].isin(columns2)] # plot only one phenotype as they can have different scales
        swarm = sns.swarmplot(data=target2, x='Variables', y='Value', hue='Taxonomy', palette=colors, ax=col_ax)
        col_ax.legend(fontsize=16) #loc='upper left')
        cm.savefig(output+'_'+pheno+'_clustermap.png')




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomic_matrix', help='Name and path of the pruned gene presence/absence file or the pruned Pfam domain presence/absence file.')
    parser.add_argument('--genomic_type', help="Data description for the title. For instance: 'Pfam domain occurence', 'gene occurence', or 'Pfam domain presence/absence'.")
    parser.add_argument('--output_prefix', help='Output prefix (can be the path).')
    parser.add_argument('--colors', help='Tab-separated color file where each line contains one of the different subspecies and the corresponding color in 3 columns (rgb).')
    parser.add_argument('--phenotype_file', help='The ; separated csv file of phenotypes.')
    parser.add_argument('--phenotype_list', help="List of phenotypes to make separate heatmaps for (separated by commas). For instance 'MIN_2000_,MIN_200_'")
    args = parser.parse_args()


    colors = load_color_scheme(args.colors)

    target = pd.read_csv(args.phenotype_file, index_col=0, sep=';')
    phenotypes = args.phenotype_list.split(',')

    genomic_matrix = pd.read_csv(args.genomic_matrix, sep=';', index_col=0)
    genomic_matrix = genomic_matrix.loc[[strain for strain in genomic_matrix.index.values if strain in target.index.values],:]


    make_dendrogram_and_clustermap(genomic_matrix, colors, args.genomic_type, args.output_prefix, target, phenotypes)

