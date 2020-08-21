import os
import numpy as np
import pandas as pd


def make_pfam_matrix(pfam_dir, out_dir, taxonomy_file, pfam_count=False, pfam_pa=False):
    """Makes and saves the dataframe(s) asked for (Pfam count-matrix and/or Pfam presence/absence-matrix)."""

    # Make sure the file name can be appended to the path later:
    if pfam_dir[-1]!='/':
        pfam_dir += '/'

    if out_dir[-1]!='/':
        out_dir += '/'

    # All the file names of PfamScan result files:
    pfam_files = os.listdir(pfam_dir)

    # Make a list of the strain names and a dictionary containing lists of their Pfams:
    strain_list = []
    pfam_strain_dict = {}
    for pfam_file in pfam_files:
        pfams = np.loadtxt(pfam_dir+pfam_file,skiprows=28, dtype='str')
        # Splitting on "." changes the the names of the NCBI RefSeq strains, preventing them from matching their IDs in the taxonomy DataFrame:
        strain_name = pfam_file.split('.')[0]
        pfam_strain_dict[strain_name] = list(pfams[:,5]) # add column 5 (Pfam IDs)
        strain_list.append(strain_name)

    print strain_list

    # Make a list of all the Pfam IDs from all the strains, where each Pfam is only listed once (for the columns of the matrix):
    combined_pfam_list = []
    for strain in pfam_strain_dict:
        combined_pfam_list += pfam_strain_dict[strain] # add all the Pfams of the strain
        combined_pfam_list = list(set(combined_pfam_list)) # a list where each Pfam is only listed once


    # Get the taxonomy of the strains:
    taxonomy = pd.read_csv(taxonomy_file, sep='; ', header=None, usecols=[0,10], engine='python')
    taxonomy.columns = ['Strains','Subspecies']
    taxonomy = taxonomy.set_index('Strains')


    # If there's no value in the last column, then there is only a ';' after the last value not a '; ' (which is what is used as sep in pd.read_csv):
    taxonomy = taxonomy.fillna('Subspecies unknown')
    taxonomy['Subspecies'] = taxonomy['Subspecies'].map(lambda x: x.rstrip(';'))


    # Create matrix or matrices:

    if pfam_count != False:
        # Make the matrix of each strain's copy number of each Pfam:
        pfam_count_matrix = pd.DataFrame(index=strain_list, columns=combined_pfam_list)
        for strain in strain_list:
            pfams_in_strain = pfam_strain_dict[strain]
            pfam_count_matrix.loc[strain] = [pfams_in_strain.count(val) for val in combined_pfam_list]
        # Add taxonomy in the first column:
        pfam_count_matrix.index.name = 'Strains'
        pfam_count_matrix = taxonomy.merge(pfam_count_matrix, on='Strains')
        # Save to file:
        pfam_count_matrix.to_csv(out_dir+pfam_count)

    if pfam_pa != False:
        # Make the matrix of each strain's presence/absence profile of the Pfams:
        pfam_pa_matrix = pd.DataFrame(index=strain_list, columns=combined_pfam_list)
        for strain in strain_list:
            pfams_in_strain = pfam_strain_dict[strain]
            pfam_pa_matrix.loc[strain] = [[0,1][val in pfams_in_strain] for val in combined_pfam_list]
        # Add taxonomy in the first column:
        pfam_pa_matrix.index.name = 'Strains'
        pfam_pa_matrix = taxonomy.merge(pfam_pa_matrix, on='Strains')
        # Save to file:
        pfam_pa_matrix.to_csv(out_dir+pfam_pa)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', help='The directory in which (only) the PfamScan results are.')
    parser.add_argument('--outdir', help='The directory to put the output files in.')
    parser.add_argument('--taxonomy_file', help='')
    parser.add_argument('--count', help='Filename to give the Pfam count matrix.')
    parser.add_argument('--pa', help='File name to give the Pfam presence/absence matrix.')
    args = parser.parse_args()

    make_pfam_matrix(args.indir, args.outdir, args.taxonomy_file, args.count, args.pa)


