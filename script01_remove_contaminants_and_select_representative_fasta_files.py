import numpy as np
import pandas as pd


def select_strains(kmerfinder_results_file,species,quast_results_file,output_file,new_files_path,fasta_file_path,min_coverage=80.0):
    if new_files_path[-1] != '/':
        new_files_path += '/'
    if fasta_file_path[-1] != '/':
        fasta_file_path += '/'

    outputlog = open(output_file+'.log','w')
    kmerfinder_results = pd.read_csv(kmerfinder_results_file, sep='\t')
    outputlog.write('Shape, all KmerFinder results:'+str(kmerfinder_results.shape)+'\n')

    with open(species,'r') as f:
        species_list = f.readline().rstrip().split('\t') #a one-line, tab-separated file containing all the allowed species.

    contaminants = set(kmerfinder_results.loc[~kmerfinder_results['Species'].isin(species_list),'Species'].values)
    if len(contaminants) > 0:
        contaminants.remove('Species')
        contaminants = np.sort(list(contaminants))
        outputlog.write('(The set of) contaminants in the samples:\n\t'+'\n\t'.join(contaminants))

    # Remove contaminants:
    kmerfinder_results = kmerfinder_results.loc[kmerfinder_results['Species'].isin(species_list)]
    outputlog.write('Shape, where species is Lactococcus lactis:'+str(kmerfinder_results.shape)+'\n')

    # Remove if coverage is too low (default is 80%):
    kmerfinder_results = kmerfinder_results.loc[pd.to_numeric(kmerfinder_results['Query_Coverage'])>float(min_coverage)]
    kmerfinder_results = kmerfinder_results.loc[pd.to_numeric(kmerfinder_results['Template_Coverage'])>float(min_coverage)]

    outputlog.write('Shape, where min coverage is above '+str(min_coverage)+'%:'+str(kmerfinder_results.shape)+'\n')
    outputlog.write('For these, the\n')
    outputlog.write('* mean Query coverage is %.2f %%\n' % (np.mean(pd.to_numeric(kmerfinder_results['Query_Coverage']))))
    outputlog.write('* mean Template coverage is %.2f %%\n' % (np.mean(pd.to_numeric(kmerfinder_results['Template_Coverage']))))


    # Select file with fewest & longest contigs based on Quast:
    quast_results = pd.read_csv(quast_results_file, sep='\t', index_col=0)
    outputlog.write("\nShape, all Quast results:"+str(quast_results.shape)+'\n')
    quast_results = quast_results.loc[quast_results.loc[:,'N50']>20000,:]
    outputlog.write("Shape, Quast results with N50 above 20000 nt:"+str(quast_results.shape)+'\n')


    kmerfinder_results = kmerfinder_results.rename(index=str, columns={'# Strain':'Assembly'})
    kmerfinder_results['Assembly'] = kmerfinder_results['Assembly'].map(lambda x: x.rstrip('.fa'))
    kmerfinder_results = kmerfinder_results.set_index('Assembly')

    sufficient_quality = kmerfinder_results.merge(quast_results, on='Assembly', validate="one_to_one")


    outputlog.write(str(len(sufficient_quality))+" files passed the quality control (species identification by KmerFinder and contig statistics (N50>20000) by Quast).\n")
    outputlog.write('Their mean Query coverage is %.2f %%' % (np.mean(pd.to_numeric(sufficient_quality.loc[:,'Query_Coverage'])))+'\n')
    outputlog.write('Their mean Template coverage is %.2f %%' % (np.mean(pd.to_numeric(sufficient_quality.loc[:,'Template_Coverage'])))+'\n')
    outputlog.write('The lowest N50 was %.1f nucleotides' % (np.min(sufficient_quality.loc[:,'N50']))+'\n')
    outputlog.write('Their mean N50 is %.1f nucleotides' % (np.mean(sufficient_quality.loc[:,'N50']))+'\n\n')
    outputlog.write('The species in the data set are:\n\t'+'\n\t'.join(set([taxonomy.split(';')[8] for taxonomy in sufficient_quality.loc[:,'Taxonomy']]))+'\n\n')
    subspecies = [taxonomy.split(';')[9] for taxonomy in sufficient_quality.loc[:,'Taxonomy']]
    outputlog.write('The subspecies are:\n\t'+'\n\t'.join([str(len([1 for sub in subspecies if sub==s]))+s for s in set(subspecies)])+'\n')
    outputlog.close()


    # Write the taxonomy of the selected strains to a separate taxonomy file:
    with open(output_file+'.taxonomy', 'w') as f:
        for file in sufficient_quality.index:
            f.write(file+'; '+sufficient_quality.loc[file,'Taxonomy']+'\n')
            new_file = open(new_files_path+file+'.fa','w')
            old_file = open(fasta_file_path+file+'.fa','r')
            for line in old_file:
                new_file.write(line)
            old_file.close()
            new_file.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--kmerfinder_results', help='A tab-separated file with the fasta file name in the first collumn, followed by the best KmerFinder match. This file must contain columns: "# Strain", "Query_Coverage", "Template_Coverage", and "Species".')
    parser.add_argument('-s','--species', help='A one-line, tab-separated file containing all the allowed species. Fasta files identified as different species will be removed.')
    parser.add_argument('-m','--min_coverage', help='The minimum coverage required of both query and template. Per default this is 80 percent.')
    parser.add_argument('-q','--quast_results', help='The transposed_report.tsv file from running Quast on all the strains.')
    parser.add_argument('-c','--output_name', help='Intermediate output file with name conversion of the selected strains. For instance "CHCC1234v4.fa,CHCC1234"')
    parser.add_argument('-f','--fasta_file_path', help='The path to the fasta files being investigated.')
    parser.add_argument('-n','--new_files_path', help='The path to the folder to put the selected (and renamed) fasta files in.')

    args = parser.parse_args()

    select_strains(args.kmerfinder_results, args.species, args.quast_results, args.output_name, args.new_files_path, args.fasta_file_path, args.min_coverage)


