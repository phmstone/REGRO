# FAQs

### What genes should I include for studying chloroplast genome degradation?
The default gene list built into the the pipeline contains 113 plastid genes, the "canonical" full set for angiosperms. 
It's best to include more genes when studying chloroplast genome degradation as this will improve the resolution at which you can see genome degradation.
Gene lists for "canonical" full sets of chloroplast genes for gymnosperms, monilophytes, mosses and hornworts, liverworts, and green algae are provided.

### Do I have to use chloroplast genome data?
No, this pipeline should work with any type of genetic data. 
Other organellar data sets, such as mitochondrial genomes, will work well.
If using nuclear data sets, it's recommended to run on and HPC for storage limitations.

### Will all the different spellings of my input gene list be recognised?
The normalisation written into the scripts handles case differences and the presence or lack of many punctuation marks.
If in doubt about alternative spellings of gene names add them to an alias file.

### I want to use the pipeline but my sequences aren't on GenBank. Will it still work?
Yes, as long as you have copies of the sequences in .fasta and .gbk format with the same name (except for the file extension). 
Put your own sequences in the folders where they would go had they been on genbank. 
It may be easier to run the pipeline through once with sequences from genbank in order to see where to put these species.

### How accurate is the updated presence/absence TSV?
The TSV output `UpdateTSV.py` from will not be an accurate assessment of whether genes found are present in full if the alignments are used straight from `aligner.py`. 
It is unlikely that a gene will be called as present if using the output from aligner.py to feed straight into `UpdateTSV.py` as it was written to be conservative and has stringent requirements in order to call a gene as "present". 
The similarity and coverage thresholds can be lowered by the user but proper checking of all the alignments is necessary to make sure that all sequences are in frame. 
The output TSV can also be edited manually.

### Can I change the text size, colours, and legend position for the heatmap figures?
Formatting changes to the heatmap figures can be made by editing the `heatMapPlot.py` script directly. 
The code chunk for colours starts and line 59, the legend position position is determined by line 112, and the code chunk for labelling starts at line 119.