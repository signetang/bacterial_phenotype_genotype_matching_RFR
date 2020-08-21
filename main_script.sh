#!/bin/bash

# NB: Certain files and directories must exist to run this script!
# It assumes the 3 directories defined bellow exist, and that ../data/_raw_data/ contains the assembled fasta files (which have long, weird names).
results="/project/results";
data="/project/data"
scripts="/project/scripts";

raw_data_path=$data"/_raw_data";
selected_file_path=$data"/selected_fasta_files";

kmerfinder_database_path="/database";
temp_KmerFinder_result=$results"/KmerFinder_result";
KmerFinder_results=$results"/KmerFinder_results.txt";

quast="/quast-5.0.2/quast.py";
quast_output_dir=$results"/quast";
allowed_species_list=$data"/species.txt";
echo "Lactococcus lactis" > $allowed_species_list;

selected_strains_overview=$results'/selected_files_based_on_KmerFinder_and_Quast';

subpopulation_level='Subspecies';
color_palette="lime,blue,deeppink";
color_scheme=$results/subspecies_colors.txt;

phenotype_file=$data/Vmax.csv;
#list of phenotype prefixes (in phenotype_file) to make separate heatmaps for:
phenotype_list='Vmax';


cd $scripts;

# Get the KmerFinder database into shared memory, to avoid reading it into memory at each run:
kma_shm -t_db $kmerfinder_database_path"/bacteria.ATG";

first_run=true; #true when running KmerFinder for the first strain in the directory

for strain in $raw_data_path/*.fa;
    do
    # Run KmerFinder on the strain:
    echo "Running KmerFinder on $strain";
    python3 /usr/src/kmerfinder.py -i $strain -o $temp_KmerFinder_result -db $kmerfinder_database_path"/bacteria.ATG" -tax $kmerfinder_database_path"/bacteria.name" -x  > /dev/null;

    if $first_run;
    then
        # Create/overwrite file for putting the best match for each strain:
        printf "# Strain\t" > $KmerFinder_results;
        # Add the header from the KmerFinder result file to the header:
        head -n 1 $temp_KmerFinder_result/results.txt >> $KmerFinder_results;
        first_run=false;
    fi;

    # Put the results for the strain into a combined file for all the strains:
    printf "%s\t" $(basename "$strain") >> $KmerFinder_results;
    head -n 2 $temp_KmerFinder_result/results.txt | tail -n 1 >> $KmerFinder_results; # The first row after the header (the best match) - or the header if no matches...

done;

# Remove temp folder (only contains the KmerFinder results from the last strain called):
rm -r $temp_KmerFinder_result/;

# Remove database from shared memory again:
kma_shm -t_db $kmerfinder_database_path"/bacteria.ATG" -destroy;

#Perform a quality control of the fasta files with Quast:
cd $raw_data_path;
python $quast -v;
python $quast --silent --output-dir $quast_output_dir --space-efficient --fast *.fa > /dev/null; #--silent doesn't seem to work, so stdout is manually sent to /dev/null
echo "Quast done.";

# Select the files which have high enough quality and are of the correct species and subspecies. Rename and put these in the folder $new_files_path:
mkdir $selected_file_path;

cd $scripts;
# For each row, the KmerFinder_results.txt file contains the name of the strain and EITHER the taxonomy or the header of the KmerFinder result.txt file.
# script02 takes this into account and only selects files, for which the Species column match the species given in the file $allowed_species_list
# (it thus dosn't match the word "Species").
python script01_remove_contaminants_and_select_representative_fasta_files.py --kmerfinder_results $KmerFinder_results --species $allowed_species_list -m 80 --quast_results $quast_output_dir"/transposed_report.tsv" --output_name $selected_strains_overview --fasta_file_path $raw_data_path --new_files_path $selected_file_path;

# Annotate strains with Prokka:
mkdir $data"/prokka";
cd $selected_file_path;

prokka --version;
for strain in *.fa;
    do
    echo "Running Prokka on "$strain;
    prokka --outdir $data"/prokka" --force --quiet --prefix $(basename $strain .fa) --locustag $(basename $strain .fa) $strain;
done;

cd $data;
# Build pangenome with Roary:
roary -f $data/roary -i 95 $data/prokka/*gff;

# Make color scheme to use for all figures in the pipeline.
#Also possible to make the tab-separated file manually: Fist column is the name of the (sub)species, the 3 next columns are the red, green, and blue values (from 0 to 1).
python $scripts/script02_color_scheme.py --color_palette $color_palette --taxonomy_file $selected_strains_overview.taxonomy --color_file $color_scheme;

# Make smaller pangenome files for each Species or Subspecies, and compute some statisitcs and make some plots for the pangenome:
# Also, make the gene presence/absence matrix file (binary: 0/1):

cd $results;
mkdir sub_pangenomes;
python $scripts/script03_pangenome_analysis.py --pangenome_file $data/roary/new_gene_presence_absence.csv --pa_matrix $data/gene_pa_matrix.txt --output_folder $results --taxonomy_file $selected_strains_overview.taxonomy --colors $color_scheme --sub_pangenome_folder $results/sub_pangenomes;

# Get the distribution of COG's in the pangenome:
python $scripts/script04_get_pangenome_COGs.py --gene_pa $data/roary/new_gene_presence_absence.csv --COG_dir $data/COG --prokka_dir $data/prokka --outfile $results/pangenome_COG_groups_whole_pangenome;

cd $results/sub_pangenomes;
for pa_file in *;
    do
    python $scripts/script04_get_pangenome_COGs.py --gene_pa $pa_file --COG_dir $data/COG --prokka_dir $data/prokka --outfile $results/pangenome_COG_groups_$pa_file;
done;

# Run PfamScan (requires the PfamScan module and an indexed Pfam HMM):
mkdir $data/Pfam;
export PERL5LIB=/PfamScan:$PERL5LIB;
export PERL5LIB=/PfamScan/Bio:$PERL5LIB;
hmmpress /Pfam27-A/Pfam-A.hmm;

cd $data/prokka;
for strain in *.faa;
    do
    echo "PfamScan on $strain ..."
    perl /PfamScan/pfam_scan.pl -outfile $data/Pfam/$strain -fasta $strain -dir /Pfam27-A;
done;

# Make Pfam matrices for input to the ML (exclude the NCBI RefSeq strains):
python $scripts/script05_make_pfam_matrix.py --indir $data/Pfam/ --outdir $data --taxonomy_file $selected_strains_overview.taxonomy --count pfam_count_matrix.txt --pa pfam_pa_matrix.txt;

# Make k-mer matrices for input to the ML (exclude the NCBI RefSeq strains)
# which exclude k-mers occuring in less than 2 % of the strains or in more than 2% of the strains:
python $scripts/script06_make_kmer_matrix.py --fasta_folder $selected_file_path --k 8 --taxonomy_file $selected_strains_overview.taxonomy --output $data/pruned_8mer_matrix.txt --count True --threshold 2
python $scripts/script06_make_kmer_matrix.py --fasta_folder $selected_file_path --k 9 --taxonomy_file $selected_strains_overview.taxonomy --output $data/pruned_9mer_matrix.txt --count True --threshold 2

# Exclude genes/Pfams which occur in less than 2 % of the strains or in more than 2% of the strains:
python $scripts/script07_prune_matrix.py --infile $data/gene_pa_matrix.txt --outfile $data/pruned_gene_pa_matrix.txt --threshold 2;
python $scripts/script07_prune_matrix.py --infile $data/pfam_count_matrix.txt --outfile $data/pruned_pfam_count_matrix.txt --threshold 2;
python $scripts/script07_prune_matrix.py --infile $data/pfam_pa_matrix.txt --outfile $data/pruned_pfam_pa_matrix.txt --threshold 2;

# Make dendrograms based on all genes/Pfams and plot heatmaps of the phenotypes in the dendrogram order:
python $scripts/script08_explore_phylogeny_as_explaining_variable.py --output_prefix $results/script08_gene_presence_absence --colors $color_scheme --genomic_type "gene presence/absence" --genomic_matrix $data/pruned_gene_pa_matrix.txt --phenotype_list $phenotype_list --phenotype_file $phenotype_file;
python $scripts/script08_explore_phylogeny_as_explaining_variable.py --output_prefix $results/script08_Pfam_presence_absence --colors $color_scheme --genomic_type "Pfam presence/absence" --genomic_matrix $data/pruned_pfam_pa_matrix.txt --phenotype_list $phenotype_list --phenotype_file $phenotype_file;

# Build and evaluate a machine learning models for the prediction of the phenotype:
python $scripts/script09_random_forest.py --genomic_matrix $data/pruned_gene_pa_matrix.txt --output_prefix $results/gene_pa --phenotype_file $phenotype_file --colors $color_scheme --phenotype_column_variables 'Measurement,Volume,Temperature,Yeast' --variable_columns ''
python $scripts/script09_random_forest.py --genomic_matrix $data/pruned_pfam_count_matrix.txt --output_prefix $results/pfam_count --phenotype_file $phenotype_file --colors $color_scheme --phenotype_column_variables 'Measurement,Volume,Temperature,Yeast' --variable_columns ''
python $scripts/script10_random_forest.py --genomic_matrix $data/pruned_8mer_matrix.txt --output_prefix $results/8mer_count --phenotype_file $phenotype_file --colors $color_scheme --phenotype_column_variables 'Measurement,Volume,Temperature,Yeast' --variable_columns ''
python $scripts/script10_random_forest.py --genomic_matrix $data/pruned_9mer_matrix.txt --output_prefix $results/9mer_count --phenotype_file $phenotype_file --colors $color_scheme --phenotype_column_variables 'Measurement,Volume,Temperature,Yeast' --variable_columns ''

# In the above four lines, the files from Supplementary Materials can be used:
# $phenotype_file = S10_File.csv
# pruned_gene_pa_matrix.txt = S6_File.txt
# pruned_pfam_count_matrix.txt = S7_File.txt
# pruned_8mer_matrix.txt = S8_File
# pruned_9mer_matrix.txt = S9_File
