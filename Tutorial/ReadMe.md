# Pipeline Tutorial

This tutorial will go through the pipeline step by step.   
Example files are included to use the pipeline to study genome degradation in 
mycoheterotrophic and autotrophic Ericaceae, and reconstruct a phylogeny.


## I. Prerequisites

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

You will also need an NCBI account in order to download sequences from genbank.
Make one [here](https://account.ncbi.nlm.nih.gov/signup/?).



## II. Downloading .gbk files from genbank

To run pipeline you will need a list of genbank IDs that relate to sequences for the taxa of interest.    
For this example we are using whole plastome sequences of some species in Ericales. The example file is called `EricalesGenbankIDs.txt`.    
Each genbank ID must be on a new line, and the list should be a `.txt` file.    
If using pipeline to investigate plastid genome degradation, make sure that all genbank IDs selected correspond to whole and not partial plastid genomes.

Use `fetchGenbank.py` to download the accession sequences into one genbank (`.gbk`) file.    

The following arguments are required:    
* `--email` The email address associated with your NCBI account.
* `--input` The name of the file containing the list of genbank IDs.
* `--output` The name of the output genbank file.    

Example command: `python3 fetchGenbank.py --email your.NCBI.account@email.address --input EricalesGenbankIDs.txt --output Ericales.gbk`   

This will create a multi-sequence genbank file in the directory `Outputs`.    



## III. Fetching and parsing annotations from genbank sequences

Now the annotations from the genbank file downloaded will be used to generate reports on which genes are present and which are absent across the taxa.    
The genes tested can be determined by the user. By default 113 genes canonically present in the angiosperm plastid genome are tested.

Use `presenceAbsence.py` to create a .TSV showing whether the "test" genes are annotated as present, absent, or pseudogenes in the downloaded genbank file.

### Inputs
Other genes or a subset of genes can be tested with the flag `gene_file`. Input a text file of gene names to bes tested, with one gene name on each line. 
An example file for ribosomal protein genes is included (`ribosomalProteinGenes.txt`).    
Most alternative spellings or formats of gene names are handled within this script's normalisation. 
If accessions you want to test contain very different names for the same gene (e.g. pafI and ycf3), then an alias file should be used.    
The alias file is a text file containing the canonical gene name and the synonym on the same line separated by a tab.    
The alias file should be called with the flag `--alias_file`. An example is included here `gene_alias.txt`.    

### Outputs
The main output is the TSV containing a row for each taxon, one column for genbank ID, one column for species name, and then a column for each gene.    
The genes are marked as present, absent, or pseudogenes according to the given annotations, where 0 = present, 1 = absent, and 2 = pseudogene.    
A text file showing presence/absence data that can be inserted into a nexus file for character state mapping in e.g. [Mesquite](https://www.mesquiteproject.org/) can be produced with `--nexus`.    
Multifasta files named by gene containing the full gene sequence for every taxon where the gene is annotated as present are generated.    
The directory containing these is named `PresentGeneMultiFastas` by default, but a different name can be chosen with the `--outdir` flag.    
Multifasta files named by gene containing the full coding sequence/gene sequence with no introns for every taxon where the gene is annotated as present are also generated.    
The directory containing these is named `PresentCodingSeqMultiFastas` by default, but a different name can be chosen with the `--outdir` flag.    

The following arguments are required:
* `--input` The .gbk file downloaded in the previous step
* `--TSV` The name of the output presence/absence matrix .TSV file.

The following arguments are optional:
* `--gene_file` A text file containing gene names separated by a new line.
* `--alias_file` A text file containing canonical gene names and synonyms separated by a tab.
* `--nexus` A text file for insertion into a nexus file for character state mapping.
* `--outdir` Output directory for whole gene sequences.
* `--coding_outdir` Output directory for exon only gene sequences.

Example command: `python3 presenceAbsence.py --input Ericales.gb --tsv EricalesPresenceAbsence.tsv --alias_file gene_alias.txt --outdir presentGeneSequences --coding_outdir --presentGeneCodingSequences`   



## IV. Visualising gene presence/absence (optional)

The output presence/absence .tsv file generated in the previous step can be hard to understand visually.    
Use `heatMapPlot.py` to create a "heatmap" style figure showing which genes are present, missing, and pseudogenised.
If not output file name is given then the image will be named     

The following argument is required:
* File name of the input .TSV with presence/absence data.

The following argument is optional:
* `-o` The name of the .png file of the plot.

Example command: `python3 heatMapPlot.py EricalesPresenceAbsence.tsv -o EricalesHeatMap.png`   


## V. 

