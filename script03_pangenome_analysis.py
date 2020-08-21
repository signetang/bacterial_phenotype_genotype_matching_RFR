import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from math import ceil
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
from script02_color_scheme import load_color_scheme

def process_pangenome(pangenome_file, pa_matrix, output_folder, taxonomy_file, colors, sub_pangenome_folder):
    """Makes a binary gene presence/absence file.
    Makes separate gene_precence_absence.csv files (gene IDs instead of 1s and 0s) for each subpopulation (species or subspecies).
    Makes histograms of number of gene groups per genome (for all strains and within each sub group).
    For each subpopulation: Finds the core, soft core, shell, cloud, and total number of genes, and the population size.
    Makes plot of the core genome and total number of genes when adding more genomes to the pangenome.

    pangenome_file should be the 'gene_presence_absence.csv' file from the Roary output.
    pa_matrix should be the file name to give the (binary) presence/absence matrix output file.
    taxonomy_file should be the file 'taxonomy_selected_files_based_on_KmerFinder_and_Quast.txt'.
    subpopulation should be either 'Subspecies' or 'Species'.
    colors should be a tab-separated file where each line contains one of the different subspecies/species and the corresponding color in 3 columns: red, green, blue.
    sub_pangenome_folder should be the folder to put the sub-pangenomes in.
    """
    if output_folder[-1] != '/':
        output_folder += '/'

    # Read taxonomy data:
    taxonomy = pd.read_csv(taxonomy_file, sep='; ', header=None, usecols=[0,10], engine='python')
    taxonomy.columns = ['Strains','Subspecies']
    taxonomy = taxonomy.set_index('Strains')
    #taxonomy.index = taxonomy.index.map(lambda x: x+'.fa') # add file extension to index names in order to merge the table with the roary pangenome matrix
    # If there's no value in the last column, then there is only a ';' after the last value not a '; ' (which is what is used as sep in pd.read_csv):
    taxonomy = taxonomy.fillna('Subspecies unknown')
    taxonomy['Subspecies'] = taxonomy['Subspecies'].map(lambda x: x.rstrip(';'))

    # Read pangenome data:
    input = pd.read_csv(pangenome_file, dtype='str', index_col=0)
    genes = input.drop(input.columns[:13].values, axis=1)
    genes = genes.T
    genes.index.rename('Strains', inplace=True)
    #genes = genes.loc[genes.index.str.contains('strain')]

    # Add the taxonomy data to the pangenome data:
    genes = taxonomy.merge(genes, on='Strains')

    # First text files:
    # Make a new pangenome file for each subpopulation:
    if sub_pangenome_folder[-1] != '/':
        sub_pangenome_folder += '/'
    subpopulation_list = list(set(genes.loc[:,'Subspecies']))
    for population in subpopulation_list:
        strain_list = genes[genes['Subspecies']==population].index
        sub_gene_matrix = input.iloc[:,:13].merge(input.loc[:,strain_list], on='Gene') #input is the original csv file
        sub_gene_matrix['No. isolates'] = np.sum(sub_gene_matrix.iloc[:,13:].fillna(0) != 0, axis=1)
        sub_gene_matrix = sub_gene_matrix.loc[sub_gene_matrix['No. isolates']!=0,:]
        sub_gene_matrix.loc[:,['No. sequences','Avg sequences per isolate','Genome Fragment','Order within Fragment','Accessory Fragment','Accessory Order with Fragment','QC','Min group size nuc','Max group size nuc','Avg group size nuc']] = ''
        sub_gene_matrix.to_csv(sub_pangenome_folder+population.split(' ')[-1]+'.csv')


    # Next text file:
    # Make a presence/absence matrix of the pangenome:
    genes_pa = genes.fillna(0)

    genes = 'empty'
    genes_pa.where(genes_pa.iloc[:,1:] == 0, 1, inplace=True)  #not the subsp. col


	# Save the presence/absence matrix (only strains for which 'strain' is part of the ID):
	#genes_pa.index = genes_pa.index.map(lambda x: x.rstrip('.fa'))
    genes_pa = genes_pa.loc[genes_pa.index.str.contains('strain')]
    genes_pa.to_csv(pa_matrix)



    # First figure:
    # Number of genes included in the pangenome after removing genes that either occur in
    # fewer strains than the threshold or are missing in fewer strains than the threshold:
    n = len(genes_pa)
    pa_counts = np.sum(genes_pa.iloc[:,1:], axis=0) # presence/absence counts per gene group
    genes_included = []
    thresholds = np.arange(51) # Makes no sense to look further than 50%. If you remove more than 50% from each end...
    for p in thresholds:
        min = n*p/100.0
        max = n*(100.0-p)/100.0

        genes_included.append(np.sum([1 if (pa_count>min and pa_count<max) else 0 for pa_count in pa_counts]))

    fig = plt.figure()
    sns.scatterplot(thresholds,genes_included,color='k')
    plt.xlim([-1,51])
    plt.ylim([int(-np.max(genes_included)/40),int(np.max(genes_included)+np.max(genes_included)/15)])
    plt.xlabel('Threshold (%)', fontsize=14)
    plt.ylabel('Number of genes included in the pangenome', fontsize=14)
    plt.title('Number of genes included in the pangenome after removing genes\nthat either occur in fewer strains than the threshold\nor are missing in fewer strains than the threshold\n', fontsize=16)
    plt.savefig(output_folder+'number_of_genes_in_pangenome_for_different_thresholds_of_occurence', bbox_inches='tight')
    plt.close()



    # Second figure:
    # Counts of gene groups (presences) in each subpopulation
    fig = plt.figure()
    ax = plt.subplot(111)
    subpopulation_list = list(set(genes_pa.loc[:,'Subspecies'])) #keeps the same order each time the set is needed
    for population in subpopulation_list:
        n = np.sum([i==population for i in genes_pa['Subspecies'].values])
        if n > 1:
            genes_pa_subpop= pd.DataFrame(data=genes_pa.loc[genes_pa['Subspecies']==population], index=genes_pa.index[genes_pa['Subspecies']==population])
            genes_pa_subpop['presence_count'] = np.sum(genes_pa.loc[genes_pa['Subspecies']==population].iloc[:,1:], axis=1)
            sns.distplot(genes_pa_subpop['presence_count'],color=colors[population],label=population, ax=ax, bins=20)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.1, box.width, box.height*0.9])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(subpopulation_list))
    plt.xlabel('Number of genes in the genome', fontsize=14)
    plt.title('Number of genes per genome', fontsize=16)
    plt.savefig(output_folder+'number_of_genes_per_genome', bbox_inches='tight')
    plt.close()


    # Next text file:
    # Pangenome stats for the entire pangenome and for each subpopulation:
    text_file = open(output_folder+'pangenome_stats.txt', 'w')
    text_file.write('\n\nPopulation\tPopulation size\tCore genes\tSoft core genes\tShell genes\tCloud genes\tTotal genes\n')

    n = len(genes_pa)
    pa_counts = np.sum(genes_pa.iloc[:,1:], axis=0) # presence/absence counts per gene group

    core = np.sum(pa_counts >= 0.99*n)
    softcore = np.sum(pa_counts >= 0.95*n) - np.sum(pa_counts >= 0.99*n)
    shell = np.sum(pa_counts >= 0.15*n) - np.sum(pa_counts >= 0.95*n)
    cloud = np.sum(pa_counts > 0) - np.sum(pa_counts >= 0.15*n)
    total = float(np.sum(pa_counts > 0))

    text_file.write('Whole pangenome\t%i\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\n' % (n, core,core/total*100, softcore,softcore/total*100, shell,shell/total*100, cloud,cloud/total*100, int(total),total/total*100))

    for population in subpopulation_list:
        n = np.sum([i==population for i in genes_pa['Subspecies'].values])
        pa_counts = np.sum(genes_pa.loc[genes_pa['Subspecies']==population].iloc[:,1:], axis=0) # presence/absence counts per gene group

        core = np.sum(pa_counts >= 0.99*n)
        softcore = np.sum(pa_counts >= 0.95*n) - np.sum(pa_counts >= 0.99*n)
        shell = np.sum(pa_counts >= 0.15*n) - np.sum(pa_counts >= 0.95*n)
        cloud = np.sum(pa_counts > 0) - np.sum(pa_counts >= 0.15*n)
        total = float(np.sum(pa_counts > 0))

        text_file.write('%s\t%i\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\t%i (%.2f%%)\n' % (population, n, core,core/total*100, softcore,softcore/total*100, shell,shell/total*100, cloud,cloud/total*100, int(total),total/total*100))

    text_file.close()




    # Third & fourth figure
    # (and their data to a text file):
    # Plot core gene count, and total gene count, as more genomes are added:
    step = 25

    number_of_strains = []
    number_of_genes = []
    number_of_core_genes = []

    n = len(genes_pa)
    for i in range(n,-1,-step):
        print 'Computing genes in core genome when', n-i+1, 'of', n, 'genomes are included...'
        include = np.array([False]*(i-1) + [True]*(n-i+1)) # include this number of strains
        for j in range(int(ceil(i/10.))):
            np.random.shuffle(include) # get a random selection of the 'this number of strains'
            number_of_strains.append(n-i+1)
            pa_counts = np.sum(genes_pa[list(include)].iloc[:,1:], axis=0) # presence/absence counts per gene group
            number_of_genes.append(np.sum(pa_counts >= 1))
            number_of_core_genes.append(np.sum(pa_counts >= ((n-i+1)*0.99) ))

    # Save info in text format:
    f = open(output_folder+'data_for_total_and_core_per_genomes_included.txt','w')
    f.write(','.join(np.array(number_of_strains).astype(str))+'\n'+','.join(np.array(number_of_genes).astype(str))+'\n'+','.join(np.array(number_of_core_genes).astype(str)))
    f.close()

    # Plot core and total gene count as more genomes are added:
    fig, ax = plt.subplots(1,1)
    sns.scatterplot(number_of_strains,number_of_core_genes, color='darkmagenta', linewidth=0, ax=ax, label='Core genes')
    sns.scatterplot(number_of_strains,number_of_genes, color='lightseagreen', linewidth=0, ax=ax, label='Total gene count')
    plt.title('Total and core genome by number of genomes included\n', fontsize=16)
    plt.ylabel('Number of genes', fontsize=14)
    plt.xlabel('Number of genomes', fontsize=14)
    plt.xlim([0,n])
    ymax = np.max(number_of_genes)
    plt.ylim([0,ymax+ymax/15.])
    plt.legend(loc=2)
    fig.savefig(output_folder+'total_and_core_per_genomes_included', bbox_inches='tight')

    # Plot (only) the core gene count as more genomes are added:
    fig, ax = plt.subplots(1,1)
    first = int(ceil(n/10.)) # number of core genomes calculated for genomes=1 (don't plot the core genomes for single genome pangenomes. That would "squish" the plot)
    sns.scatterplot(number_of_strains[first:],number_of_core_genes[first:], color='k', linewidth=0, ax=ax, label='Whole pangenome')
    xmax = n # max value on the x-axis

    # And the core genomes of the sub-populations:
    for population in subpopulation_list:
        # overwrite for each subspecies:
        number_of_strains = []
        number_of_core_genes = []
        n = np.sum([i==population for i in genes_pa['Subspecies'].values])
        if n>step:
            for i in range(n-step,-1,-step):
                print 'Computing genes in core genome of',population,'when', n-i+1, 'of', n, 'genomes are included...'
                include = np.array([False]*(i-1) + [True]*(n-i+1)) # include this number of strains
                for j in range(int(ceil(i/10.))):
                    np.random.shuffle(include) # get a random selection of the 'this number of strains'
                    pa_counts = np.sum(genes_pa.loc[genes_pa['Subspecies']==population][list(include)].iloc[:,1:], axis=0) # presence/absence counts per gene group
                    #pa_counts = np.sum(genes_pa[list(include)].iloc[:,1:], axis=0) # presence/absence counts per gene group
                    number_of_core_genes.append(np.sum(pa_counts >= ((n-i+1)*0.99) ))
                    number_of_strains.append(n-i+1)
            sns.scatterplot(number_of_strains, number_of_core_genes, color=colors[population], linewidth=0, ax=ax, label=population)

    plt.title('Core genes by number of genomes included', fontsize=16)
    plt.ylabel('Number of core genes', fontsize=14)
    plt.xlabel('Number of genomes', fontsize=14)
    plt.xlim([0,xmax])
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+box.height*0.15, box.width, box.height*0.85])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=1)
    fig.savefig(output_folder+'core_by_genomes_included', bbox_inches='tight')





if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pangenome_file', help="The 'gene_presence_absence.csv' file from the Roary output.")
    parser.add_argument('--pa_matrix', help='Path and name to give the (binary) presence/absence matrix output file.')
    parser.add_argument('--output_folder', help='Path for the output files (plots and text files).')
    parser.add_argument('--taxonomy_file', help="Name and path of the file 'taxonomy_selected_files_based_on_KmerFinder_and_Quast.txt'")
    parser.add_argument('--colors', help='Tab-separated file where each line contains one of the different subspecies/species and the corresponding color in 3 columns: red, green, blue.')
    parser.add_argument('--sub_pangenome_folder', help='(Existing) folder to put the sub-pangenomes in.')
    args = parser.parse_args()


    colors = load_color_scheme(args.colors)
    process_pangenome(args.pangenome_file, args.pa_matrix, args.output_folder, args.taxonomy_file, colors, args.sub_pangenome_folder)
