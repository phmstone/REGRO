# 🌿🧬🌷 REGRO 🌷🧬🌿

## <ins>R</ins>etrieval <ins>E</ins>ngine for <ins>G</ins>ene <ins>R</ins>ecovery in <ins>O</ins>rganelles

REGRO is a suite of Python scripts for identifying and extracting genes and gene fragments (pseudogenes) from organellar genomes.
Using an input list of GenBank IDs, the annotations will be parsed to produce alignments of present genes and a report of gene presence/absence/pseudogenisation (GenBank mode).
Further scripts are included to find genes and gene fragments that may have been missed by the annotations on the GenBank files (discovery mode).


Go to the [tutorial](https://github.com/phmstone/REGRO/tree/main/Tutorial) for a guide on how to run REGRO.
If you have questions, please check the [FAQ page](https://github.com/phmstone/REGRO/blob/main/FAQs.md) before raising an issue.


REGRO is designed for users who are comfortable in the command line environment.
Some useful tutorials for biologists getting started using the command line can be found [here](https://github.com/mossmatters/introToCmdLine/blob/master/introToCmdLine.pdf) and [here](https://datacarpentry.github.io/shell-genomics/). 
Work through these first if you feel you need more experience before using REGRO.

### Pipeline diagram

![Flowchart diagram for REGRO showing inputs, outputs, script order, and optional steps](/Images/pipeline.png)


### Supported platforms 

REGRO will has been tested on macOS and Linux. Windows users can use REGRO with Windows Subsystem for Linux (WSL).


### Installation 

You can obtain REGRO in one of the following ways.

**Clone the repository** (recommended)    
```
git clone https://github.com/phmstone/REGRO.git    
cd REGRO
```

**Download as a ZIP**    
Click the green Code button on the GitHub repository page.
Select Download ZIP.
Extract the archive and navigate to the directory.

**Command line download**    
```
wget https://github.com/phmstone/REGRO/archive/refs/heads/main.zip    
unzip main.zip
```

### Prerequisites 

**Required dependencies**  
* BLAST
* MAFFT 

These are only needed to run REGRO in discovery mode.    
Both can be installed with [Conda](https://anaconda.org/), [Homebrew](https://brew.sh/), or manually and then added to your $PATH. 
If you are planning to run REGRO on an HPC then check if these tools are available as modules, as installation may be unnecessary.    

**Required third party Python packages**  
* Bio
* certifi
* numpy
* matplotlib
* pandas
* seaborn

Python packages can be installed with [pip](https://pypi.org/project/pip/), as shown in the [tutorial](https://github.com/phmstone/REGRO/tree/main/Tutorial).

**Using Conda to download dependencies**  

If you prefer to use Conda then a yaml file is also provided in the home directory.

```
name: REGRO

channels:
  - conda-forge
  - bioconda

dependencies:
  - python=3.12
  - biopython=1.87
  - blast=2.17.0
  - mafft=7.525
  - certifi=2026.7.22
  - matplotlib=3.11.1
  - numpy=2.5.1
  - pandas=3.0.5
  - seaborn=0.13.2
```

This will ensure no version issues are encountered with the different dependencies needed for REGRO.

Use it like this:

```
git clone https://github.com/username/TOOL.git
cd REGRO
conda env create -f environment.yml
conda activate REGRO
```
