# FAQs

### What genes should I include for studying genome degradation?
The default gene list built into the the pipeline contains 113 plastid genes, the "canonical" full set for angiosperms. 
It's best to include more genes when studying genome degradation as this will improve the resolution at which you can see pseudogenisation and gene loss.
Gene lists for "canonical" full sets of genes for monilophytes, bryophytes, and green algae plastid genomes and the angiosperm mitochondrial genome are provided.

### Do I have to use chloroplast genome data?
No, this pipeline should work with any type of genetic data. 
Other organellar data sets, such as mitochondrial genomes, should work well.
If using nuclear data sets, it's recommended to run on an HPC for storage reasons. 
If your reads are not yet assembled, try using [PhyloHerb](https://github.com/lmcai/PhyloHerb) instead.

### Will all the different spellings of my input gene list be recognised?
The normalisation written into the scripts handles case differences and the presence or lack of many punctuation marks.
If in doubt about alternative spellings of gene names, add them to an alias file.

### How many alternative spellings/names can be added to the alias file?
As many as you like.

### What if my "reference" sequence doesn't include all the genes I want to study?
Choose a different sequence that does contain the full complement of genes, or remove genes missing from the reference from the gene list.

### I want to use the pipeline but my sequences aren't on GenBank. Will it still work?
Yes, as long as you have copies of the sequences in .fasta and .gbk format with the same name (except for the file extension). 
Put your own sequences in the folders where they would go had they been on genbank. 
It may be easier to run the pipeline through once with sequences from genbank in order to see where to put your own.

### How accurate is the updated presence/absence TSV?
The TSV output `UpdateTSV.py` from will not be an accurate assessment of whether genes found are present in full if the alignments are used straight from `aligner.py`. 
It is unlikely that a gene will be called as present if using the output from aligner.py to feed straight into `UpdateTSV.py` as it was written to be conservative and has stringent requirements in order to call a gene as "present". 
The similarity and coverage thresholds can be lowered by the user but proper checking of all the alignments is necessary to make sure that all sequences are in frame. 
The output TSV can also be edited manually to reflect changes in gene status discovered by manual editing.

## Does this test if detected tRNA genes have functional anticodons?
No. The pipeline uses similarity and percentage coverage of the reference sequence only to call non-coding protein genes' status. 
I recommend you look into using [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) in addition to pipeline if you need to know whether your the tRNA genes have functional anticodons.

### Can I change the text size, colours, and legend position for the heat map figures?
Formatting changes to the heat map figures can be made by editing the `heatMapPlot.py` script directly. 
The code chunk for colours starts and line 59, the legend position position is determined by line 112, and the code chunk for labelling starts at line 119.