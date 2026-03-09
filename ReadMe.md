# 🌿🧬🌷 Plastid genome content project 🌷🧬🌿

This pipeline is a suite of python scripts designed to investigate plastid genome degradation. 
Using an input list of GenBank IDs, the annotations will be parsed to produce alignments of present genes and a report of gene presence/absence/pseudogenisation.
Further scripts are included to find genes and gene fragments that may have been missed by the annotations on the GenBank files.
Go to the [tutorial](https://github.com/phmstone/PlastidGenomeContent/tree/main/Tutorial) for a guide on how to run it.
If you have questions, please check the [FAQ page](https://github.com/phmstone/PlastidGenomeContent/blob/main/Tutorial/FAQs.md) before raising an issue.


The pipeline is designed for users who are comfortable in the command line environment.
Some useful tutorials for biologists getting started using the command line can be found [here](https://github.com/mossmatters/introToCmdLine/blob/master/introToCmdLine.pdf) and [here](https://datacarpentry.github.io/shell-genomics/). 
Work through these first if you feel you need more experience before using this pipeline.


## Installation

You can obtain the pipeline in one of the following ways.

**Clone the repository** (recommended)    
```
git clone https://github.com/phmstone/PlastidGenomeContent.git    
cd PlastidGenomeContent
```

**Download as a ZIP**    
Click the Code button on the GitHub repository page.
Select Download ZIP.
Extract the archive and navigate to the directory.

**Command line download**    
```
wget https://github.com/phmstone/PlastidGenomeContent/archive/refs/heads/main.zip    
unzip main.zip
```

## Necessary python modules for the pipeline to run
* SeqIO (from Bio)
* Seq (from Bio.Seq)
* SeqRecord (from Bio)
* argparse
* os
* time
* subprocess
* certifi
* Entrez (from Bio)
* numpy
* Path (from pathlib)
* re
* collections
* seaborn
* certifi
* matplotlib
* pandas

Python packages can be installed with [pip](https://pypi.org/project/pip/), as shown in the [tutorial](https://github.com/phmstone/PlastidGenomeContent/tree/main/Tutorial).

## Software necessary for full pipeline
* BLAST
* MAFFT
* Exonerate    

All can be installed with [Conda](https://anaconda.org/), [Homebrew](https://brew.sh/), or manually and then added to your $PATH. 
If you are planning to run this on an HPC then check if these tools are available as modules, as installation may be unnecessary.
If installing Exonerate with Homebrew, use the bioinformatics specific channel before installing.

## Inputs needed
* List of "test" GenBank IDs in .txt file, one ID per line
* List of reference GenBank IDs that contain all plastid genes in .txt file, one ID per line

The rest of the files needed will be generated innately.






## Potential errors (Graham lab can ignore!)

`Failed for NC_068711: <urlopen error [SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self signed certificate in certificate chain (_ssl.c:992)>`
Had to put "import certifi" and "os.environ["SSL_CERT_FILE"] = certifi.where()" for this to get fixed.