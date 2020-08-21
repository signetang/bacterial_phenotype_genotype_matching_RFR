# KmerFinder database:
#version=latest; #not the latest due to errors...
version=20190108_stable/bacteria;
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/$version/bacteria* -P database;
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/$version/config -P database;
tar -xzf database/bacteria.tar.gz -C database;
rm database/*.tar.gz;

# Pfam-A Hidden Markov Model:
mkdir Pfam27-A;
cd Pfam27-A;
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz;
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.dat.gz;
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/active_site.dat.gz;
gunzip Pfam-A.hmm.dat.gz;
gunzip Pfam-A.hmm.gz;
gunzip active_site.dat.gz;
cd ..

mkdir COG;
cd COG;
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab;
wget ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab;
