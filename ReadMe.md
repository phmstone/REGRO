# 🌿🧬🌷 REGRO 🌷🧬🌿

## <u>R</u>etrieval <u>E</u>ngine for <u>G</u>ene <u>R</u>ecovery in <u>O</u>rganelles

REGRO is a suite of python scripts designed to investigate plastid genome degradation. 
Using an input list of GenBank IDs, the annotations will be parsed to produce alignments of present genes and a report of gene presence/absence/pseudogenisation (GenBank mode).
Further scripts are included to find genes and gene fragments that may have been missed by the annotations on the GenBank files (discovery mode).
Go to the [tutorial](https://github.com/phmstone/PlastidGenomeContent/tree/main/Tutorial) for a guide on how to run REGRO.
If you have questions, please check the [FAQ page](https://github.com/phmstone/PlastidGenomeContent/blob/main/Tutorial/FAQs.md) before raising an issue.


REGRO is designed for users who are comfortable in the command line environment.
Some useful tutorials for biologists getting started using the command line can be found [here](https://github.com/mossmatters/introToCmdLine/blob/master/introToCmdLine.pdf) and [here](https://datacarpentry.github.io/shell-genomics/). 
Work through these first if you feel you need more experience before using REGRO.


### Installation

You can obtain REGRO in one of the following ways.

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

### Required third party python packages
* Bio
* certifi
* numpy
* matplotlib
* pandas
* seaborn

Python packages can be installed with [pip](https://pypi.org/project/pip/), as shown in the [tutorial](https://github.com/phmstone/PlastidGenomeContent/tree/main/Tutorial).

### Required software
* BLAST
* MAFFT 

These are only needed to run REGRO in discovery mode.    
Both can be installed with [Conda](https://anaconda.org/), [Homebrew](https://brew.sh/), or manually and then added to your $PATH. 
If you are planning to run REGRO on an HPC then check if these tools are available as modules, as installation may be unnecessary.    



