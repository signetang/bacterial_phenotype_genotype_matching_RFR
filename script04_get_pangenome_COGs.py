import os
import numpy as np
import pandas as pd
from math import ceil

def COG_annotations_of_pangenome(roary_pangenome_presence_absence_file, COG_dict, min_included=0, max_included='all'):

	df = pd.read_csv(roary_pangenome_presence_absence_file, dtype='str')
	# Do not include the RefSeq strains:
	strains = df.columns[14+35:] # the first cols are: "Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"

	if max_included == 'all':
		max_included = len(strains)

	pangenome = df.loc[[(min_included<=int(i) and int(i)<=max_included) for i in df.loc[:,'No. isolates'].values],:]

	#pangenome = df.iloc[:,14:].fillna(0)
	#pangenome = df.loc[[(min_included<=int(i) and int(i)<=max_included) for i in np.sum(pangenome != 0,axis=1)],:]

	COG_categories = []

	# Get list of the COG categories of the gene groups:
	for row in range(len(pangenome)): #for each gene group
		genes = list(set(pangenome.iloc[row,14:])) # list of unique genes in the gene group
		COGs = []
		# For each gene in the gene group, get the COG annotation (if any):
		for gene in genes:
			if gene in COG_dict:
				COGs.append(COG_dict[gene])

		#if len(set(COGs)) > 1: # if a group has more than one COG annotation
		#	print 'There was more than one COG assigned to the gene group:', pangenome.iloc[row,0], set(COGs)

		if len(COGs) > 0:
			COG_categories.append(max(set(COGs), key=COGs.count)) # Append the gene group's (most frequently occuring) COG

		#for category in list(set(COGs)):
		#	COG_categories.append(category) # Include all the gene group's COGs, not only the most frequent

	# Count occurences of each category:
	COG_categories = [i+': '+str(COG_categories.count(i)) for i in set(COG_categories)]

	COG_categories = ['This part of the pangenome contains '+str(len(pangenome))+' gene groups.\nThose annotated with COGs are distributed in following categories:']+COG_categories

	for category in COG_categories:
		print category

	return COG_categories


def make_COG_dict(prokka_dir, COG_ID_to_func_dict):
	""" Returns a dictionary of gene IDs to COGs for all genes which prokka annotated with COGs. """

	print '\nPlease wait. Computing dictionary of COGs...'

	if prokka_dir[-1] != '/':
		prokka_dir += '/'

	# Get list of .tsv files from the Prokka output:
	prokka_files = os.listdir(prokka_dir)
	prokka_tsv_files = []
	for file in prokka_files:
		if file[-4:] == '.tsv':
			prokka_tsv_files.append(file)

	# Make dictionary of gene IDs to COG category:
	COG_dict = {}

	for tsv_file in prokka_tsv_files:
		df = pd.read_csv(prokka_dir+tsv_file, sep='\t', usecols=[0,5])
		df = df.dropna() # remove rows without COG ID
		for i in range(len(df)):
			COG_dict[df.iloc[i,0]] = COG_ID_to_func_dict[df.iloc[i,1]] # Get the functional category of the COG.
	print 'Done computing the dictionary.'
	return COG_dict


def make_COG_ID_to_func_dict(COG_dir):
	"""Returnes dictionary from COG id to functional category (for instance from 'COG0623' to 'Lipid transport and metabolism')."""
	if COG_dir[-1] != '/':
		COG_dir += '/'

	categories = pd.read_csv(COG_dir+'fun2003-2014.tab', sep='\t')
	catdict = {}
	for i in range(len(categories)):
		catdict[categories.iloc[i,0]] = categories.iloc[i,1]

	functions = pd.read_csv(COG_dir+'cognames2003-2014.tab', sep='\t')
	dict = {}
	for i in range(len(functions)):
		dict[functions.iloc[i,0]] = catdict[functions.iloc[i,1][0]] # [0] as there can be more than one category, and the first is considered the main one.
	return dict


def analyse_pangenome(gene_pa, COG_id_to_category_dict, output_file):
	"""Finds which COG categories the pangenome has in it's core, soft core, shell, and cloud genomes, and in the entire pangenome.
	Prints the output and stores it in a file (output_file).
	Takes three inputs:
	gene_pa: The (path and) name of the gene_presence_absence.csv file from the Roary pangenome
	COG_id_to_category_dict: A dictionary from COG IDs to categories
	output_file: The (path and) name of the output file to store the results in."""

	output_file = output_file.rstrip('csv').rstrip('.')

	# Get the number of strains:
	df = pd.read_csv(gene_pa, dtype='str')
        strains = df.columns[14:]
	n_strains = len(strains)

	#print '\nStrict core genome (genes in 100% of strains):'
	#COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict, n_strains, n_strains)

	# Find out how the COG categories distribute over the pangenome's core, soft core, shell, and cloud genomes, as well as the entire pangenome:
	print '----------------------------------------------'
	print 'COG group distribution in the Prokka pangenome'
	print '----------------------------------------------'
	print '\nCore genome (genes in >99% of strains):'
	core = COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict, ceil(0.99*n_strains), n_strains)
	print '\nSoft core genome (genes in >95% and <99% of strains):'
	soft_core = COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict, ceil(0.95*n_strains), ceil(0.99*n_strains)-1)
	print '\nShell genome (genes in >15% and <95% of strains):'
	shell = COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict, ceil(0.15*n_strains), ceil(0.95*n_strains)-1)
	print '\nCloud genome (genes in <15% of strains):'
	cloud = COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict, 1, ceil(0.15*n_strains)-1)
	print '\nEntire pangenome (all genes in the pangenome):'
	all = COG_annotations_of_pangenome(gene_pa, COG_id_to_category_dict)

	f = open(output_file+'.txt','w')
	f.write('Core:\n'+"\n".join(core)+'\n\nSoftcore:\n'+"\n".join(soft_core)+'\n\nShell:\n'+"\n".join(shell)+'\n\nCloud:\n'+"\n".join(cloud)+'\n\nEntire pangenome:\n'+"\n".join(all))
	f.close()



if __name__ == "__main__":
        import argparse
        parser = argparse.ArgumentParser()
	parser.add_argument('--gene_pa', help='The gene_presence_absence.csv file from the Roary output.')
        parser.add_argument('--COG_dir', help='The directory containing the two files "fun2003-2014.tab" and "cognames2003-2014.tab".')
        parser.add_argument('--prokka_dir', help='The Prokka output directory. The script uses the .tsv files.')
	parser.add_argument('--outfile', help='Path and filename of file to write the output to.')
        args = parser.parse_args()

	# Make dict translating COG ID to functional category:
	COG_func_dict = make_COG_ID_to_func_dict(args.COG_dir)
	# Make dict translating gene ID to COG functional category:
	COG_id_to_category_dict = make_COG_dict(args.prokka_dir, COG_func_dict)
	# Find which COG categories the pangenome has in it's core, soft core, shell, and cloud genomes, and in the entire pangenome,
	# and write the output to a file (outfile):
	analyse_pangenome(args.gene_pa, COG_id_to_category_dict, args.outfile)

