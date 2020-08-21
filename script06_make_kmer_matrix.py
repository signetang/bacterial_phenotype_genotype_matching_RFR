#import numpy as np
import os
import pandas as pd
import numpy as np
from script07_prune_matrix import prune_genes


def reverse_DNA(sequence):
    """Returns the reverse-complement sequence.
    Only use this function if it has been tested, that the sequences only contain the codes: A,C,G,T,N"""
    sequence = sequence.replace('T','U')
    sequence = sequence.replace('A','T')
    sequence = sequence.replace('U','A')
    sequence = sequence.replace('G','U')
    sequence = sequence.replace('C','G')
    sequence = sequence.replace('U','C')
    return sequence[::-1]


def kmer_counter(sequence, k, kmers, count):
    """Counts the k-mers of a string ('sequence') and adds the counts to the dictionary ('kmers').
    Or, if 'count' is set to False, the function reports presence/absence of the k-mers."""

    for start_pos in range(len(sequence)-k+1):
        kmer = sequence[start_pos:start_pos+k]
        reverse = reverse_DNA(kmer)  ### Get the opposite strand k-mer
        kmer = np.sort([kmer, reverse])[0] #### Choose the aplphabetically first version
        if kmer not in kmers:
            kmers[kmer] = 1
        elif count:
            kmers[kmer] += 1
    return kmers


def contig_reader(fasta_file, k, kmers, count):
    """Counts the k-mers of all contigs in a fasta file. Adds them to the dictionary 'kmers'.
    Or, if 'count' is set to False, adds the presence/absence of k-mers in all the contigs to the dictionary 'kmers'."""

    contig = False
    with open(fasta_file,'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                if contig != False: # if it is not the first contig
                    # Add the k-mers of the previously read contig:
                    kmers = kmer_counter(contig, k, kmers, count)
                ####new_contig = True
                # Initiate new contig:
                contig = ''
            else:
                # Add sequence to contig:
                contig += line.rstrip()
        # Add the k-mers of the last contig:
        kmers = kmer_counter(contig, k, kmers, count)

    return kmers


def kmer_matrix(folder, k, taxonomy_file, output, count, threshold):
    """Counts the k-mers of each fasta file in a folder. Makes a matrix of the counts or of the presence/absence if 'count' is set to False.
    The matrix will have a row for each strain. The first column is for the taxonomy and the following columns are for each k-mer."""

    if folder[-1] != '/':
        folder += '/'

    kmer_matrix = {}
    kmerlist = []
    i = 1
    for fasta_file in os.listdir(folder):
        # Only include strains whith names starting by 'strain':
        if fasta_file[:6] == 'strain':
            kmers = {}
            # Count all the k-mers in the fasta file:
            kmers = contig_reader(folder+fasta_file, k, kmers, count)
            # Add new k-mers to the list:
            kmerlist += [s for s in kmers]
            kmerlist = list(set(kmerlist))
            print i, len(kmerlist)
            i += 1
            kmer_matrix[fasta_file] = kmers

    kmer_matrix = pd.DataFrame(kmer_matrix)
    kmer_matrix = kmer_matrix.fillna(0)
    kmer_matrix = kmer_matrix.astype(int).transpose()
    kmer_matrix.index = kmer_matrix.index.map(lambda x: x.rstrip('.fa'))


    # Add the taxonomy of the strains:

    # Read taxonomy data:
    taxonomy = pd.read_csv(taxonomy_file, sep='; ', header=None, usecols=[0,10], engine='python')
    taxonomy.columns = ['Strains','Subspecies']
    taxonomy = taxonomy.set_index('Strains')
    # If there's no value in the last column, then there is only a ';' after the last value not a '; ' (which is what is used as sep in pd.read_csv):
    taxonomy = taxonomy.fillna('Subspecies unknown')
    taxonomy['Subspecies'] = taxonomy['Subspecies'].map(lambda x: x.rstrip(';'))

    # Add taxonomy in the first column:
    kmer_matrix.index.name = 'Strains'
    kmer_matrix = taxonomy.merge(kmer_matrix, on='Strains')

    # Prune the matrix:
    kmer_matrix = prune_genes(kmer_matrix, threshold, threshold)

    # Save to file:
    kmer_matrix.to_csv(output, sep=';')


def nucleotides_acceptable(folder):
    """Tests if the fasta files contain any other nucleic acid codes than A, C, G, T, or N.
    Returns True if these are the only codes."""

    acceptable = ['A','C','G','T','N']

    strains = [strain for strain in os.listdir(folder) if strain[:6]=='strain']

    for strain in strains:
        letters = ''
        with open('../data/selected_fasta_files/'+strain,'r') as fil:
            for lines in fil:
                if lines[0] != '>':
                    line = lines.rstrip()
                    letters += line
        for l in letters:
            if l not in acceptable:
                return False

    return True



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_folder', help='The path and folder of the fasta files (containing only the fasta files).')
    parser.add_argument('--k', help='The length of the k-mers.')
    parser.add_argument('--taxonomy_file', help='The path and name of the KmerFinder result file.')
    parser.add_argument('--output', help='The new file to write the k-mer matrix to.')
    parser.add_argument('--count', help="'True' to return the counts of k-mers or 'False' to return the presence/absence.")
    parser.add_argument('--threshold', help='Percentage. Genes either occuring or missing in fewer than this percentage of strains are not included in the pruned matrix.')
    args = parser.parse_args()


    if nucleotides_acceptable(args.fasta_folder): # test that the sequences only contain the codes: A,C,G,T,N.

        if args.count in ['True','true','T','t']:
            count = True
        else:
            count = False

        print 'Making matrix of k-mer counts of the files in '+args.fasta_folder
        kmer_matrix(folder=args.fasta_folder, k=int(args.k), taxonomy_file=args.taxonomy_file, output=args.output, count=count, threshold=int(args.threshold))

        print 'Done. Stored in '+args.output

    else:
        print '\nFailed test.\nOnly use this script if the sequences only contain the codes: A,C,G,T,N.\n'
