import pandas as pd
from math import floor,ceil


def presences(gene_profile):
    return sum([1 for i in gene_profile if i>0])


def prune_genes(genes, more_than=10, less_than=10):
    """Removes columns from the genes pd.DataFrame if the sum is >more_than and <less_than.
    Identical columns are collapsed so that only one instance of the identical columns are kept.
    The names of the identical columns are concatteneted."""

    # Recalculate from percentages to number of strains:
    n_strains = len(genes.index)
    more_than = int(floor(more_than/100.*n_strains))
    less_than = int(ceil((1.0-less_than/100.)*n_strains))

    print 'Non pruned shape:', genes.values.shape
    print 'Genes (or Pfams) occurring in less than',more_than,'strains or in more than',less_than,'strains are excluded from the pruned matrix.'
    #print '\nGenes (or Pfams) with identical occurrence profiles across all the strains:'

    # Save taxonomy information in first column of the pruned matrix:
    pruned = pd.DataFrame()
    pruned['Taxonomy'] = genes.iloc[:,0]

    # Drop the taxonomy column and convert the remaining genes dataframe to integers:
    genes = genes.iloc[:,1:]
    genes = genes.astype(int)

    i = 0

    for gene in list(genes):
        # The ocurence profile of either genes or Pfams:
        profile = genes.loc[:,gene] # The 1's and 0's (or copy number)

        # Adds the domain if it occurs the right amount of times
        pres = presences(profile)
        if pres>more_than and pres<less_than:

            # Subtracts the profile from all previously added profiles, gets the sum of the absolute values of the profiles.
            abs_diff = pruned.iloc[:,1:].sub(profile, axis='index').abs().sum(axis=0)

            # If the sum of a profile is 1, then the new profile is identical to one previously added:
            if list(abs_diff ).count(0) > 0:

                # If the new profile is identical to a previously added one, concatenate them:
                previous_gene = abs_diff[abs_diff==0].index[0]
                #print genes.loc[:,[previous_gene,gene]],'\n'
                #print previous_gene+','+gene
                pruned = pruned.rename(index=str,columns={previous_gene:previous_gene+','+gene})

            else: # Else, add the new, distinct profile:
                pruned[gene] = profile

       # Remove column from original matrix to save space??:
       genes = genes.drop(gene, axis=1)
       i += 1

    print 'Pruned shape:', pruned.values.shape, '\n'

    return pruned


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', help='Comma-separated table (with header) of occurences of genes or Pfams in the strains. Rows are for different strains, collumns for different genes.')
    parser.add_argument('--outfile', help='Name to give the resulting file containing the pruned matrix.')
    parser.add_argument('--threshold', help='Percentage. Genes either occuring or missing in fewer than this percentage of strains are not included in the pruned matrix.')
    args = parser.parse_args()

    print 'Pruning matrix: '+args.infile
    input_table =  pd.read_csv(args.infile, sep=',', index_col=0)
    pruned = prune_genes(input_table, int(args.threshold), int(args.threshold))
    pruned.to_csv(args.outfile, sep=';')



