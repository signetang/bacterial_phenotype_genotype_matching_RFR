import os
import numpy as np
import pandas as pd
from collections import Counter


def make_geneid2group(pangenome_file):
    """ Makes a dictionary of gene_id to gene_group."""
    pangenome = pd.read_csv(pangenome_file, dtype='str', index_col=0)
    #pangenome = pangenome.fillna('')

    n_rows, n_cols = pangenome.shape
    gene_groups = pangenome.index.values

    geneid2descriptions = {}
    for row in range(n_rows):
        geneid2descriptions[gene_groups[row]] = gene_groups[row]+', '+', '.join([str(n) for n in pangenome.iloc[row,[0,1]]])

    geneid2group = {}
    for row in range(n_rows):
        for col in range(13+35,n_cols,1):
            gene_id = pangenome.iloc[row,col]
            if pd.notna(gene_id):
                if '\t' in gene_id:
                    for paralog in gene_id.split('\t'):
                        geneid2group[paralog] = gene_groups[row]
                else:
                    geneid2group[gene_id] = gene_groups[row]



    return geneid2group, geneid2descriptions


def reverse_DNA(sequence):
    """Get the reverse transcript."""
    sequence = sequence.replace('T','U')
    sequence = sequence.replace('A','T')
    sequence = sequence.replace('U','A')
    sequence = sequence.replace('G','U')
    sequence = sequence.replace('C','G')
    sequence = sequence.replace('U','C')
    return sequence[::-1]


def strand_dict(gene_file):
    strand_df = pd.read_csv(gene_file[:-3]+'gff', sep='\t', error_bad_lines=False, comment='#',header=None)
    strand_df = strand_df.iloc[:,[6,8]]
    strand_df = strand_df.dropna()
    strand_df.columns = ['strand','ID']
    strand_df.ID = strand_df.ID.map(lambda x: x.split(';')[0].split('=')[1])
    strand_df = strand_df.set_index('ID')

    strand = {}
    for i in range(len(strand_df)):
        strand[strand_df.index[i]] = strand_df.iloc[i,0]

    return strand


def get_gene_from_kmer(kmer,gene_file):
    """ Makes a list of those genes in the gene_file which contains the kmer.
    Goes through the lines of the nucleotide gene fasta files from Prokka.
    Makes each gene into a single string. Checks if the kmer is in the string.
    If so, adds the gene ID to the list which is returned. """

    strand_of = strand_dict(gene_file)
    gene = 'initial'
    strand_of['initial'] = 'none'

    with open(gene_file, 'r') as file:
        gene_list = []
        sequence = '' # initial sequence which does not contain the kmer

        for line in file:
            # When reading a fasta header, update the gene name and initiate new sequence:
            if line[0] == '>':
                # If previous gene contained the kmer, add the gene id to the list:
                if kmer in sequence: # look for the k-mer if the gene is on the + strand
                    gene_list.append(gene)
                elif reverse_DNA(kmer) in sequence: # look for the reverse k-mer if the gene is on the - strand
                    gene_list.append(gene)
                # Initiate new gene:
                gene = line.split(' ')[0][1:]
                sequence = ''
            else:
                sequence += line.rstrip()

    return gene_list


def get_gene_groups_from_kmers(kmer_list, gene_file_dir, geneid2group, geneid2descriptions):
    """Makes a dicitonary of the gene groups which each kmer in the kmer_list occur in."""

    if gene_file_dir[-1] != '/':
        gene_file_dir += '/'

    # Make a dictionary where the kmers are keys to lists (for the gene_groups)
    gene_groups_of_kmer = {}
    for kmer in kmer_list:
        gene_groups_of_kmer[kmer] = []

    gene_file_list = [file for file in os.listdir(gene_file_dir) if (file.endswith('.ffn') and file.startswith('strain'))]

    for kmer in kmer_list:
        for gene_file in gene_file_list:
            gene_list = get_gene_from_kmer(kmer, gene_file_dir+gene_file)
            if len(gene_list) > 0:
                gene_groups = [geneid2group[gene] for gene in gene_list if gene in geneid2group] # It seems that roary threw out some genes from the pangenome.. Hence the 'if'
                gene_groups_of_kmer[kmer] += gene_groups

    return gene_groups_of_kmer


def get_gene_groups_from_pfams(pfam_list, pfam_file_dir, geneid2group, geneid2descriptions):
    """Prints the names of gene groups in the pangenome, which have a gene that contains the given Pfam domain."""

    if pfam_file_dir[-1] != '/':
        pfam_file_dir += '/'

    # Make a dictionary where the pfams are keys to lists (for the gene_groups)
    gene_groups_of_pfam = {}
    for pfam in pfam_list:
        gene_groups_of_pfam[pfam] = []

    pfam_file_list = [file for file in os.listdir(pfam_file_dir) if file.startswith('strain')]

    for pfam_file in pfam_file_list:
        pfam2geneid = pd.read_csv(pfam_file_dir+pfam_file, delim_whitespace=True, usecols=[0,5], index_col=1, skiprows=27, header=None) # the pfam col becomes the index col

        for pfam in pfam_list:
            # If the Pfam is found in the organism, add the corresponding gene group name(s) to the list
            if pfam in pfam2geneid.index.values:
                for gene in pfam2geneid.loc[pfam,].values.flatten():
                    if gene in geneid2group:
                        gene_groups_of_pfam[pfam].append(geneid2group[gene])
                    else:
                        print gene,

    return gene_groups_of_pfam


def get_sequence(pangenome_file, gene_groups, fasta_folder, fasta_suffix='.faa', n_sequences=1):
    pangenome = pd.read_csv(pangenome_file, dtype='str', index_col=0)

    if fasta_folder[-1] != '/':
        fasta_folder += '/'

    genomes = pangenome.columns
    fasta = ''
    for gene_group in gene_groups:
        counter = 1

        if n_sequences == 'all':
            n_sequences = len(genomes) # this may be higher than the total number of sequences
        for i in range(13,len(genomes)): #go through all columns
            if pd.notna(pangenome.loc[gene_group,genomes[i]]): #when encountering a gene
                for gene_ID in pangenome.loc[gene_group,genomes[i]].split('\t'):
                    fasta += get_from_fasta(fasta_folder+genomes[i]+fasta_suffix, gene_ID, new_name=gene_group+'__'+gene_ID) #get the sequence of the gene
                    counter += 1
            if counter > n_sequences:
                break # break loop as soon as enough sequences have been extracted (but always get all in a paralog group)

    return fasta


def get_from_fasta(file, gene, new_name=''):
    sequence = ''
    record = False
    with open(file,'r') as f:
        for line in f:
            if line[0] == '>':
                if line[1:len(gene)+1] == gene:
                    record = True
                    if len(new_name) > 0:
                        sequence += '>'+new_name+'\n'
                    else:
                        sequence += line
                elif record:
                    # If record is True when encountering a new sequence (line[0] == '>'), then the entire protein has been read. Return it:
                    return sequence
            elif record:
                # If the correct gene fasta header has been encountered, then add the list to the sequence:
                sequence += line
         # If at the end of the file:
        if record:
            return sequence
        else:
            print new_name, gene
            return 0

    raise Exception('No sequence found for '+gene+' in '+file)


def get_contig_from_kmer(kmer,genome_dir):
    """ Makes a list of those contigs in each genome_file which contains the kmer.
    Goes through the lines of the nucleotide genome fasta files from Prokka.
    Makes each contig into a single string. Checks if the kmer is in the string.
    If so, adds the contig ID and position to a strain dictionary which is returned. """
    if genome_dir[-1] != '/':
        genome_dir += '/'

    genome_list = [file for file in os.listdir(genome_dir) if (file.endswith('.fna') and file.startswith('strain'))]

    positions = {}

    for strain in genome_list:
        strain_name = strain.split('/')[-1].split('.')[0]

        with open(genome_dir+strain, 'r') as file:
            contig_list = []
            sequence = '' # initial sequence which does not contain the kmer

            for line in file:
                # When reading a fasta header, update the contig name and initiate new sequence:
                if line[0] == '>':
                    # If previous gene contained the kmer, add the gene id to the list:
                    if kmer in sequence:
                        start = sequence.index(kmer)
                        stop = sequence.index(kmer)+len(kmer)
                        contig_list.append([start,stop,contig,sequence[start-5:stop+5]])
                        # Initiate new gene:
                    contig = line.rstrip()[1:]
                    sequence = ''
                else:
                    sequence += line.rstrip()
        if contig_list != []:
            positions[strain_name] = contig_list

    return positions


print 'Started...'

# Make a dictionary of gene_id to gene_group:
geneid2group, geneid2descriptions = make_geneid2group('../../data/roary/new_gene_presence_absence.csv')
print 'Pangenome dictionary done..'

# Print Pfam name and the counts of genes from gene group in which the Pfam is found:
pfam_list = ['PF00082.17', 'PF02225.17', 'PF06280.7', 'PF00746.16']
gene_groups_of_pfam = get_gene_groups_from_pfams(pfam_list, '../../data/Pfam/', geneid2group, geneid2descriptions)
for pfam in pfam_list:
    print '\n', pfam
    print pd.Series([geneid2descriptions[gene_group] for gene_group in gene_groups_of_pfam[pfam]]).value_counts()


# Print Pfam name and the counts of genes from gene group in which the Pfam is found:
pfam_list = ['PF00639.16','PF02502.13','PF02254.13','PF02386.11','PF13493.1','PF02669.10','PF03814.10','PF02702.12','PF07274.7','PF11361.3','PF01432.15','PF07083.6','PF12730.2','PF05649.8','PF02486.14','PF04103.10','PF04967.7','PF02796.10']
gene_groups_of_pfam = get_gene_groups_from_pfams(pfam_list, '../../data/Pfam/', geneid2group, geneid2descriptions)
for pfam in pfam_list:
    print '\n', pfam
    print pd.Series([geneid2descriptions[gene_group] for gene_group in gene_groups_of_pfam[pfam]]).value_counts()


# 9-mers Combined + and - strand:
kmer_list = ['AGGGCCCAG', 'CCGAGACCG','AAGATCTAC', 'GCCGTCGAC', 'AGAGTCCGG', 'CCGGGTAGC', 'GCCAGGGAC', 'CTACCCGGC', 'CGCGGCGTA', 'CGCCGCGGC', 'CGCAGCCCC', 'ACTGGGCCC', 'GATCTAACC', 'CGGAACCCG', 'CGACCTACA', 'AGCCCTATG', 'CACGCTGCG']
gene_groups_of_kmer = get_gene_groups_from_kmers(kmer_list, gene_file_dir='../../data/prokka/', geneid2group=geneid2group, geneid2descriptions=geneid2descriptions)
for kmer in kmer_list:
    print '\n', kmer, '/',  reverse_DNA(kmer)
    print pd.Series([geneid2descriptions[gene_group] for gene_group in gene_groups_of_kmer[kmer]]).value_counts()


# 8-mers Combined + and - strand:
kmer_list = ['AACCGGGG', 'CCTGGCCA', 'CGTATACG', 'CCGGGTAG', 'CCCGGCCC', 'GCCGGGTA', 'CCCGCGGG', 'AGACCGGG', 'CGCGCCCC', 'CCGATCGA', 'CGGCCCCG', 'CGGCCCGA', 'GCCGTGCC', 'CGGTCCGA', 'GCGCGCCC', 'CCCACGGG', 'GGGGCCCC']
gene_groups_of_kmer = get_gene_groups_from_kmers(kmer_list, gene_file_dir='../../data/prokka/', geneid2group=geneid2group, geneid2descriptions=geneid2descriptions)
for kmer in kmer_list:
    print '\n', kmer, '/',  reverse_DNA(kmer)
    print pd.Series([geneid2descriptions[gene_group] for gene_group in gene_groups_of_kmer[kmer]]).value_counts()

