# Pipeline Tutorial

This tutorial will go through the pipeline step by step.   
Example files are included to use the pipeline to study genome degradation in 
mycoheterotrophic and autotrophic Ericaceae, and reconstruct a phylogeny.


## I. Installation and prerequisites


You can obtain the pipeline in one of the following ways.

Clone the repository (recommended)    
```
git clone https://github.com/phmstone/PlastidGenomeContent.git    
cd REPOSITORY
```

Download as a ZIP     
Click the Code button on the GitHub repository page.
Select Download ZIP.
Extract the archive and navigate to the directory.

Command line download       
```
wget https://github.com/phmstone/PlastidGenomeContent/archive/refs/heads/main.zip    
unzip main.zip
```

Install the following software to run the pipeline:
* [Python 3](https://www.python.org/downloads/)
* [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html)
* [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)

If running on an HPC (e.g. Digital Research Alliance of Canada) then pre-installed modules can be used.     
Be sure to load them before attempting to run Pipeline.

Several third party Python packages are needed:
* biopython
* numpy
* certifi
* pandas
* matplotlib
* seaborn

Run `pip install pandas matplotlib seaborn biopython numpy certifi` to install these if necessary.

You will also need an NCBI account in order to download sequences from GenBank.
Make one [here](https://account.ncbi.nlm.nih.gov/signup/?).



## II. Downloading .gb files from GenBank

To run pipeline you will need a list of GenBank IDs that relate to sequences for the taxa of interest.    
For this example we are using whole plastome sequences of some species in Ericales. The example file is called `EricalesGenBankIDs.txt`.    
Each GenBank ID must be on a new line, and the list should be a `.txt` file.    
If using pipeline to investigate plastid genome degradation, make sure that all GenBank IDs selected correspond to whole and not partial plastid genomes.

Use `fetchGenBank.py` to download the accession sequences into one GenBank (`.gb`) file.    

The following arguments are required:    
* `--email` The email address associated with your NCBI account.
* `--input` The name of the file containing the list of GenBank IDs.
* `--output` The name of the output GenBank file.    

**Example command:** 
```
python3 fetchGenBank.py --email your.NCBI.account@email.address --input EricalesGenBankIDs.txt --output Ericales.gb
```   

This will create a multi-sequence GenBank file in the directory `Outputs`.    



## III. Fetching and parsing annotations from GenBank sequences

Now the annotations from the multi-sequence GenBank file downloaded in step II will be used to generate reports on which genes are present and which are absent across the taxa.    
The genes tested can be determined by the user. By default 113 genes canonically present in the angiosperm plastid genome are tested.

Use `presenceAbsence.py` to create a .TSV showing whether the "test" genes are annotated as present, absent, or pseudogenes in the downloaded GenBank file.

#### Inputs
Other genes or a subset of genes can be tested with the flag `gene_file`. Input a text file of gene names to bes tested, with one gene name on each line. 
An example file for ribosomal protein genes is included (`ribosomalProteinGenes.txt`).    
Most alternative spellings or formats of gene names are handled within this script's normalisation. 
If accessions you want to test contain very different names for the same gene (e.g. _pafI_ and _ycf3_), then an alias file should be used.    
The alias file is a text file containing the canonical gene name and the synonym on the same line separated by a tab.    
The alias file should be called with the flag `--alias_file`. An example alias file is included [here](https://github.com/phmstone/PlastidGenomeContent/blob/main/Tutorial/exampleFiles/gene_alias.txt).    

#### Outputs
The main output is the TSV file containing a row for each taxon, one column for GenBank ID, one column for species name, and then a column for each gene tested.    
The genes are marked as present, absent, or pseudogenes according to the given annotations from the input GenBank file, where 0 = present, 1 = absent, and 2 = pseudogene.    
A text file showing presence/absence data that can be inserted into a nexus file for character state mapping in e.g. [Mesquite](https://www.mesquiteproject.org/) can be produced with `--nexus`.    
Multifasta files named by gene containing the full gene sequence for every taxon where the gene is annotated as present are generated.    
The directory containing these is named `PresentGeneMultiFastas` by default, but a different name can be chosen with the `--outdir` flag.    
Multifasta files named by gene containing the full coding sequence/gene sequence with no introns for every taxon where the gene is annotated as present are also generated.    
The directory containing these is named `PresentCodingSeqMultiFastas` by default, but a different name can be chosen with the `--outdir` flag.    

The following arguments are required:
* `--input` The .gb file downloaded in the previous step
* `--tsv` The name of the output presence/absence matrix .TSV file.

The following arguments are optional:
* `--gene_file` A text file containing gene names separated by a new line.
* `--alias_file` A text file containing canonical gene names and synonyms separated by a tab.
* `--nexus` A text file for insertion into a nexus file for character state mapping.
* `--outdir` Output directory for whole gene sequences. This directory is produced automatically, but the name can be changed with this flag.    
* `--coding_outdir` Output directory for exon only gene sequences.
* `--coding_outdir` Output directory for sequences annotated as pseudogenes.    

**Example command:** 
```
python3 presenceAbsence.py --input Ericales.gb --tsv EricalesPresenceAbsence.tsv --alias_file gene_alias.txt --outdir presentGeneSequences --coding_outdir presentGeneCodingSequences
```   



## IV. Visualising gene presence/absence (optional)

The output presence/absence .tsv file generated in the previous step can be hard to understand visually.    
Use `heatMapPlot.py` to create a "heatmap" style figure showing which genes are present, missing, and pseudogenised.
If no output file name is given then the image will be named `heatMapPlot.png`.     

The following argument is required:
* `--input` File name of the input .TSV with presence/absence data.

The following argument is optional:
* `--output` The name of the .png file of the plot.

**Example command:** 
```
python3 heatMapPlot.py --input EricalesPresenceAbsence.tsv --output EricalesHeatMap.png
```

## V. Finding other sequences with BLAST

Sometimes not all the genes present in a genome are successfully annotated, this is especially true for pseudogenes. With BLAST we can uncover genes that may have been missed by the original annotation method, or not included in the final version of the genome on GenBank for some other reason.     
The script `blastPresenceAbsence.py` does this.
To save time and computing power, this script only uses BLAST to try and pull out genes that are annotated as pseudogenes or are missing on GenBank. Genes annotated as present are assumed to be annotated correctly.    
This script downloads fasta files for the plastid genomes of all sequences with at least one gene annotated as missing/pseudogenised, and makes databases for BLAST from these fasta files.    
BLAST is performed using gene sequences from "reference sequences" as queries.

#### Inputs
Reference sequences should be supplied to act as query sequences for the BLAST searches. The reference sequences should have functional copies of all the genes you are testing. Ideally they are closely related taxonomically.    
An example file containing some GenBank IDs for sequences containing all 113 of the default angiosperm plastid genes is included (`referenceIDs.txt`)    
The reference sequences should be identified by their GenBank IDs and the reference sequence file should have one GenBank ID per line in a .txt file.    
If some of the sequences supplied in the original GenBank ID list for testing contain functional copies of all the genes, then this script will enable them to be used as references for BLAST.

#### Outputs
Fasta files from the GenBank sequences in the original .txt file list, BLAST directories, BLAST results, and .gb files for the "reference" sequences are all produced and organised automatically in directories.     
Fasta files, BLAST database files, and .gb files are all named after the relevant GenBank ID. BLAST results are stored in directories named after the GenBank ID that was used to make the BLAST database, and individual results files are named by gene name.    
By default, a new directory called `Blast` is created for all the output files and directories to go into, but the directory name can be changed.

The following arguments are required:
* `--input` The presence/absence .TSV file downloaded in the previous step
* `--email` The email address associated with your NCBI account.

The following arguments are optional:
* `--reference-ids` A text file containing GenBank IDs known to have all genes of interest separated by a new line. Optional because if one of the IDs in the input .TSV contains all genes of interest, that GenBank ID can be used as a "reference".
* `--blast-type` Choose between tblastx or blastn (defualt is blastn).
* `--outdir` A directory that all ouputs are put into (default name is `Blast/ReferenceGeneSequences`).    

**Example command:** 
```
python3 blastPresenceAbsence.py --input EricalesPresenceAbsence.tsv --email your.NCBI.account@email.address --reference-ids referenceIDs.txt --blast-type blastn
```   


## VI. Parsing BLAST results

The results from the BLAST searches can be used to extract putative gene sequences, which can then be put into alignments with the present genes that were extracted in `presenceAbsence.py` to determine if these sequences are full genes or gene fragments/pseudogenes.    
Pseudogenes annotated on GenBank are not extracted and put into the alignments here, as they should have sufficient homology to be pulled out with this approach.   
Sequences with very long hits where the entire SSC region of a plastid may have been captured are flagged and put in a file called `FilesToCheckAgain.txt`.
In order to include all sequences found as "present" for your taxa of interest, use `--present-genes` and point to the directory containing those multifastas.
Gene/taxon combinations where no BLAST hit was found are written out in a file called `NoHits.txt`.
There are two methods available for parsing BLAST results `blastProcessing.py` or `blastProcessing-singleSeq.py`. 
`blastProcessing.py` merges hits based on their absolute start and end coordinates in the genome, where as `blastProcessing-singleSeq.py` merges hits based on the genomic distance between the hits and will only append the single longest merged sequence.   
Read the [FAQ](https://github.com/phmstone/PlastidGenomeContent/blob/main/Tutorial/FAQs.md) if you're not sure which one to use. 


`blastProcessing.py`
The following arguments are required:
* `--output-dir` A directory to store unaligned multifastas of hits from the Blast results.
* `--blast-dir` Directory containing BLAST results.
* `--reference-dir` Directory containing reference gene multifasta files.
* `--genome-dir` Directory containing whole plastid fasta sequences. 

The following arguments are optional:
* `--flanking-region` The number of base pairs taken before and after the the start and end of the blast hit. the default value is `0`.
* `--ir-cutoff` On rare occasions the script will incorrectly process multiple blast hits together to create a very long sequence. Files with a sequence over the length set here will be flagged. The defualt value is `5000`.
* `--present-genes` Directory containing "present" sequences (full gene or coding sequence) for genes of interest in fasta format. 

**Example command:** 
```
python3 blastProcessing.py --output-dir unalignedMultiFastas --ir-cutoff 4000 --blast-dir Blast/Results  --reference-dir Blast/ReferenceGeneSequences --genome-dir Blast/PlastidSequences
``` 


`blastProcessing-singleSeq.py`
The following arguments are required:
* `--output-dir` A directory to store unaligned multifastas of hits from the Blast results.
* `--blast-dir` Directory containing BLAST results.
* `--reference-dir` Directory containing reference gene multifasta files.
* `--genome-dir` Directory containing whole plastid fasta sequences. 

The following arguments are optional:
* `--merge-gap` The maximum number of base pairs in between to blast hits for those hits to be merged into one sequence.
* `--present-genes` Directory containing "present" sequences (full gene or coding sequence) for genes of interest in fasta format. 

**Example command:** 
```
python3 blastProcessing-singleSeq.py --output-dir unalignedMultiFastas --merge-gap 800 --blast-dir Blast/Results  --reference-dir Blast/ReferenceGeneSequences --genome-dir Blast/PlastidSequences
``` 


## VII. Aligning

In order to determine if the sequences extracted are pseudogenes or full genes likely to be functional copies, the multifastas generated in the previous step are aligned.
Multifastas are aligned with MAFFT to produce a variety of alignments with `aligner.py`. 
Nucleotide, codon, and protein alignments, as well as an unaligned protein sequence multifasta file are output for all protein coding genes.
Alignments containing only the reference sequences are also produced for each gene.
Gene names starting with "trn" or "rrn" are assumed to be non protein coding, so only nucleotide alignments are produced for these genes.

The following arguments are required:
* `--input` The directory containing unaligned multifasta files produced in the previous step.
* `--output` The directory where the aligned multifasta files will be stored.

**Example command:** 
```
python3 aligner.py --input unalignedMultiFastas --output alignedMultiFastas
``` 



## VIII. Re-evaluating gene presence/absence

`updateTSV.py` will assess alignments to determine whether the "found" sequences for each taxon-gene combination represent a pseudogene or full gene copy.
Any blast hit put into an alignment will count as a pseudogene. 
To be called as a present gene, the sequence must meet minimum coverage and similarity thresholds (adjustable) to the reference sequences, and not contain any premature stop codons.
For genes starting with "trn" or "rrn" the minimum coverage and similarity thresholds alone are used.
A change log file is output, where each line shows the old and new values and reason for change for each gene-taxon combination changed in the presence/absence .TSV.

The following arguments are required:
* `--ogTSV` The original presence/absence .TSV created from GenBank annotations.
* `--alignDir` The directory containing aligned multifasta files produced in the previous step.

The following arguments are optional:
* `--outTSV` The updated presence/absence TSV file name (default is `updatedPA.tsv`)
* `--changeLog` The file name for the log TSV recording changes made (default is `changeLog.tsv`)
* `--minCov` The minimum coverage needed to a reference sequence for a gene to be called as present (default is `0.95`).
* `--minSim` The minimum similarity needed to a reference sequence for a gene to be called as present (default is `0.95`).

**Example command:** 
```
python3 updateTSV.py --ogTSV EricalesPresenceAbsence.tsv --alignDir alignedMultiFastas
```



## IV. Visualising the updated gene presence/absence (optional)

Use `heatMapPlot.py` again to create a "heatmap" style figure showing which genes are present, missing, and pseudogenised with the updated presence/absence .tsv.
Make sure that the figure has a different name to the one created in step IV with the GenBank annotations.

The following argument is required:
* `--input` File name of the updated .TSV with presence/absence data.

The following argument is optional:
* `--output` The name of the .png file of the plot.

**Example command:** 
```
python3 heatMapPlot.py --input updatedPA.tsv --output EricalesHeatMap-updated.png
``` 



## Pipeline Complete!



