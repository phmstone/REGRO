# FAQs

### What genes should I include for studying genome degradation?
The default gene list built into the the pipeline contains 113 plastid genes, the "canonical" full set for angiosperms. 
It's best to include more genes when studying genome degradation as this will improve the resolution at which you can see pseudogenisation and gene loss.
Gene lists for "canonical" full sets of genes for monilophytes, bryophytes, and green algae plastid genomes and the angiosperm mitochondrial genome are provided.

### Do I have to use chloroplast genome data?
No, this pipeline should work with any type of genetic data. 
Other organellar data sets, such as mitochondrial genomes, should work well provided that the whole genome is assembled as one sequence with one GenBank ID.
This could work for nuclear datasets, it's recommended to run on an HPC for storage reasons.
If your reads are not yet assembled, try using [PhyloHerb](https://github.com/lmcai/PhyloHerb) instead.

### Will all the different spellings of my input gene list be recognised?
The normalisation written into the scripts handles case differences and the presence or lack of many punctuation marks.
If in doubt about alternative spellings of gene names, add them to an alias file.

### How many alternative spellings/names can be added to the alias file?
As many as you like.

### How will I know if I need to use an alias file?
If one of the genomes you are looking into was annotated a while ago, it may have alternative names for some genes like _ycf1_ or _ycf2_. 
If you've already made a presence/absence table or figure, then see if whole suites of genes are missing unexpectedly (e.g. tRNA, rRNA, _ndh_). 
This could indicate a systemic formatting difference in these gene names that the normalisation is unable to handle.

### How do I know if the sequences I want to use as references have all genes in my gene list annotated as present?
Use `fetchGenbank.py` and `presenceAbsence.py` with an input of the gene list you want to use and a list of GenBank IDs for the taxa you are considering using as references. 
The outputs of `presenceAbsence.py` will show which taxa contain all genes in your gene list and are suitable to be used as references.

### What if my "reference" sequence doesn't include all the genes I want to study?
Choose a different sequence that does contain the full complement of genes, or remove genes missing from the reference from the gene list.    
You can use the first two scripts with a large number of GenBank IDs as input to identify accessions that do have all genes. 
The easiest way to see this by making a heat map figure, or by summing up the total count of the presence/absence .tsv by accession. Accessions with all genes present will have a sum of 0.

### One of taxa I am studying has a much lower number of present genes than I am expecting, what could be wrong?
It's possible that this GenBank accession is not very well annotated, so even though genes are present they are not annotated and will not be recognised by the first part of the pipeline. 
Alternative gene names may also have been used to annotate this accession, have a look at the GenBank file and make an alias file if this is the case.

### I want to use the pipeline but my sequences aren't on GenBank. Will it still work?
Yes, as long as you have copies of the sequences in .fasta and .gbk format with the same name (except for the file extension). 
Put your own sequences in the folders where they would go had they been on genbank. 
It may be easier to run the pipeline through once with sequences from genbank in order to see where to put your own.

### Will this pipeline recover genes with trans-spliced introns?
In angiosperm chloroplast genomes, _rps12_ is the only gene with trans-spliced introns. 
Alignments for genes with trans-spliced will take more manual curation, but if `blastProcessing.py` is used then the full sequence should be recovered.
If `blastProcessing-singleSeq.py` is used, then the whole gene will not be recovered, as this script extracts the longest sequence for each taxon the different exon(s) could be extracted from different sequences. 
If the different exons for the gene are located in different chromosomes (i.e. have different GenBank IDs) then this script will skip over those sequences and not recover the genes. 
This should only be a problem for mitochondrial sequences.

### How accurate is the updated presence/absence TSV?
The TSV output `UpdateTSV.py` from will not be an accurate assessment of whether genes found are present in full if the alignments are used straight from `aligner.py`, manual inspection and editing is likely necessary. 
It is unlikely that a gene will be called as present if using the output from aligner.py to feed straight into `UpdateTSV.py` as it was written to be conservative and has stringent requirements in order to call a gene as "present". 
The similarity and coverage thresholds can be lowered by the user but proper checking of all the alignments is necessary to make sure that all sequences are in frame. 
The output TSV can also be edited manually to reflect changes in gene status discovered by manual editing. 
Take extra care interpreting these results if you are working with taxa that use an [alternative genetic code](https://www.pnas.org/doi/10.1073/pnas.1816822116) or have [high rates of RNA editing](https://academic.oup.com/nar/article/31/9/2417/1080299?login=false).    

### Does this test if detected tRNA genes have functional anticodons?
No. The pipeline uses similarity and percentage coverage of the reference sequence only to call non-coding protein genes' status. 
I recommend you look into using [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE) in addition to pipeline if you need to know whether your the tRNA genes have functional anticodons.

### Can I change the text size, colours, and legend position for the heat map figures?
Formatting changes to the heat map figures can be made by editing the `heatMapPlot.py` script directly. 
The code chunk for colours starts and line 59, the legend position position is determined by line 112, and the code chunk for labelling starts at line 119.