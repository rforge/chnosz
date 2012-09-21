# the following data files support calculations using the 
# RefSeq database (release 49, 2011-09-07)
protein_refseq.csv: overall (average) amino acid composition of all proteins for each
  microbial genome in the RefSeq collection (n=3471)
taxid_names.csv: taxid, phylum name and species name for 3471 microbial taxa

# these functions/scripts have the following purpose (output files listed in parentheses):
mkfaa.sh - combine gzipped sequence files into one big FASTA file (refseq49.faa)
gencat.sh - extract gi number, taxid, sequence length from RefSeq release catalog (gi.taxid.txt)
protein.refseq.R - get average amino acid composition for each taxid from gzipped sequence files (protein_refseq.csv)
taxid.names.R - get taxonomic names for each taxid represented (taxid_names.csv)

# bash scripts assume a GNU/Linux-like operating system
# timings were made for processing RefSeq 49 on a recent (2009) intel laptop

# get the list of files and entries in the database
1. download 'release49.files.installed' and 'RefSeq-release49.catalog.gz' from NCBI
   (ftp://ftp.ncbi.nih.gov/refseq/release/release-catalog)
2. gzip -d RefSeq-release49.catalog.gz [1.7 GB]

# download stuff
3. list URLS for the microbial protein sequence files:
     grep microbial.*.protein.faa* release49.files.installed | \
       sed -e "s/^/ftp\:\/\/ftp.ncbi.nih.gov\/refseq\/release\/microbial\//g" > urllist
4. download the files using 'wget -i urllist' [1340 files, 2.0 GB]
5. move the .gz files to a directory named 'protein'
6. run ls protein/*.gz > filelist
7. use 'mkfaa.sh' to combine the sequences into a single file 'refseq49.faa' [4.0 GB, ~3 minutes]

# protein stuff
8. use 'gencat.sh' to generate gi.taxid.txt from RefSeq-release49.catalog
   note that the intermediate file gi.taxid.unsrt may have to be edited manually 
     -- see instructions in gencat.sh
   when done, the output of 'cat gi.taxid.txt | wc -l'  
   should be equal to 'grep "^>" refseq49.faa | wc -l'
   (for microbial proteins in RefSeq 49, the number is 10368471)
9. generate protein_refseq.csv in R:  [~5.5 hours]
   > source("protein.refseq.R")
   > protein.refseq()
   note that this depends on gi.taxid.txt and the .faa.gz files in the 'protein' directory

# taxonomy stuff
10. edit 'taxid.names.R' so that 'taxdir' points to the directory where the files
    'names.dmp' and 'nodes.dmp' are present. these files can be downloaded from
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz (accessed on 2011-09-06)
11. source 'taxid.names.R' to generate the file 'taxid_names.csv' [~1.5 hours]
