# Scripts

## fetchGenbank.py
This script takes an input of a list of genbank IDs in a text file and makes an output file of one big genbank file. At the moment that big genbank file goes into the `Outputs` directory within `Scripts`.  
The file name of the output file has to be given as an argument on the command line. An email address associated with an account on genbank/NCBI must also be given to run the script.  

## presenceAbsence.py
This script converts the annotations present in genbank files to a TSV with a column for each gene and a row for each taxon.  
The value entered represents the status of the gene in that organism (0 = present, 1 = missing, 2 = pseudogene).  

## presentGeneMultiFasta.py
This script makes multi fastas by gene for all chloroplast genes (could be manually changed for desired genes of interest) where the gene is annotated as being present in the genbank file.  
It creates a directory called `PresentGeneMultiFastas` where these multifasta files are kept, but the name of this directory can be chosen using `--outdir` on the command line.  

## blastPresenceAbsence.py
This script uses the .TSV generated from `presentGeneMultiFasta.py` and a list genbank IDs to use as plastid gene reference sequences input by the user.  
The presence/absence profiles from the TSV are read in by gene for each species. Species from this TSV that have all genes present will be added to the reference species list.  
The reference species list should ideally contain genbank IDs of whole plastid sequences of species closely to the taxa of interest that contain all plastid genes of interest.  
For each species in the TSV, blast searches are performed on genes that are categorised as either missing or pseudogenised. The blast searches are done with a reference sequence from the user's input reference species list.  
The script outputs a directory `Blast` with four subdirectories; `Databases`, `PlastidSequences`, `ReferenceGeneSequences`, `ReferenceGenomes`, and `Results`.  
* `Databases` database files for every genbank accession with missing genes or pseudogenes that are needed for blast
* `PlastidSequences` fasta files for every genbank accession with missing or pseudogenes
* `ReferenceGeneSequences` A fasta file for every gene of interest created from the first reference sequence
* `ReferenceGenomes` Genbank files for all reference sequences (including those added after checking presence/absence profiles)
* `Results` Text files from the output of each blast search performed to look for missing or pseudo genes using a reference sequence

## blastProcessing.py
This script uses the blast results generated from `blastPresenceAbsence.py` to generate multifasta files for all of the genes tested.  
The multifasta files for each gene contain sequences from the reference plastid sequences as outlined in `blastPresenceAbsence.py` and sequences from the test plastid genomes that were identified as having high similarity to the gene with blast searches.  
The sequences in the multifasta file from the test plastid genomes are labelled with the sequence coordinates and the reference sequence that was used to find that specific hit.  
Test plastid - gene combinations that had no hits found with blast are output to `NoHitsFiles.txt`.  
Exceptionally long blast hits are are not written out to the multifasta files, but the .txt blast result file is flagged and written to `filesToCheckAgain.txt`.  
Mandatory inputs are directories containing blast results, reference gene sequences, and test plastid fasta files.   
All of the above directories are generated inside `Blast` by `blastPresenceAbsence.py` automatically.  
An output directory name can be specified, or the multifasta files and .txt flagging files will be put in a directory called `unalignedMultifastas`.  
The threshold for how long a blast hit has to be to be flagged can be adjusted with `--ir-cutoff` but the default is 5000.  
The number of base pairs flanking the blast hit to be included in the final multifasta file can be specified hit with `--flanking-region`, but the default is set to 0.  

Things to do 
* fix flagging/ strangeness around IR in `blastProcessing.py`
* make a list of all modules needed for pipeline
* deal with genes that have introns
* make an alignment script
* make a script that checks if genes are pseudogenes (completeness, indels, premature stop codon)

