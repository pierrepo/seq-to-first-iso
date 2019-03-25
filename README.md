# Seq-to-first-iso

> Compute first isotopologues intensity from peptide sequence

## Contents
- [Background](#background)
- [Installation](#installation)
  - [With a virtual environment](#with-a-virtual-environment)
  - [Manually](#manually)
- [Usage](#usage)
  - [Options](#options)
  - [Examples](#examples)
- [Credits](#credits)

## Background

[(Back to top)](#contents)

<!--GET IMAGE OF SPECTRO ? One of the main challenge of mass spectrometry is the identification of peptides
Isotopologues are molecules that differ only in their isotopic composition  
**seq-to-first-iso** aims to provide a way to compute M0 and M1 in 12C conditions with unlabelled amino acids  
More details here or other section ?
This project is based on... (SLIM)
SAY THAT ABUNDANCES ARE C ONLY, AND WHAT THEY ARE (establish convention: M0_NC, 12C, formula_X)-->

## Installation

[(Back to top)](#contents)

**Note**: this is the installation process for LINUX machines, depending on your operating system, the commands might differ

### With a virtual environment

This tutorial uses [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "conda installation guide") to set up the packages in a virtual environment

- Download _environment.yml_
- After going in the directory containing _environment.yml_, create the environment with:

```shell
$ conda env create --file environment.yml
```
You should now have a conda environment named SLIM
To go into the environment, type:

```shell
$ conda activate SLIM
```
You should now be able to use the script inside the environment

### Manually

**Download the following prerequesites:**

- [Python](https://www.python.org/downloads/) 3.7 or higher
- [pyteomics](https://pyteomics.readthedocs.io/en/latest/) 4.0 or higher
- [pandas](https://pandas.pydata.org/) 0.24 or higher


## Usage

[(Back to top)](#contents)

**If you installed with conda, make sure to activate the conda environment**  
The script takes a file with one sequence of amino acids per line and returns a tsv of the file with columns:

|sequence|mass|formula|formula_X| M0_NC | M1_NC | M0_12C | M1_12C |
|--------|----|-------|---------|-------|-------|--------|--------|

The script can be called with:

```shell
$ python seq-to-first-iso.py filename [-o output_name] [-n amino_acids...]
```
Optional arguments are in square brackets  
This will create _filename.tsv_ if filename is a correct file


### Options

- `-h, --help`:  
Provide a help page

- `-o, --output`:  
Change the name of the output file

- `-n, --non_labelled_aa`:  
Take 1 or more amino acid separated by a comma

### Examples

- You can provide a list of amino acids which will keep default isotopic abundance:

Supposing _peptides.txt_ :

```
YAQEISR
VGFPVLSVKEHK
LAMVIIKEFVDDLK
```

The command
```shell
$ python seq-to-first-iso.py peptides.txt -n V,W
```
will create _peptides.tsv_ :

|sequence| mass| formula| M0_NC| M1_NC| M0_12C| M1_12C|
|--------|-----|--------|------|------|-------|-------|
YAQEISR| 865.42938099921| C37H59O13N11| 0.6206414140575179|	0.280870823368276| 0.9206561231798033| 0.05161907174495234|
VGFPVLSVKEHK| 1338.7659712609| C63H102O16N16| 0.4550358985377136| 0.34506032928190855| 0.7589558393662944| 0.18515489894512063|
LAMVIIKEFVDDLK| 1632.91606619252| C76H128O21N16S1| 0.36994021481230627| 0.3373188347614264| 0.7475090558698947| 0.15292723586285323|

Where, in 12C enrichment conditions, the isotopologue intensity M0_12C and M1_12C are computed with unlabelled Valine and Tryptophan  (V and W have default isotopic abundance)


- You can change the name of the output file:

```shell
$ python seq-to-first-iso.py peptides.txt -o sequence
```
will create a file named _sequence.tsv_

## Credits

[(Back to top)](#contents)

- **Bioconda**:
  - Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. “Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences”. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

- **pyteomics**:
  - Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

  - Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717

- **MIDAs**:
  - Alves G, Ogurtsov AY, Yu YK (2014) Molecular Isotopic Distribution Analysis (MIDAs) with adjustable mass accuracy. J Am Soc Mass Spectrom, 25: 57-70. DOI: 10.1007/s13361-013-0733-7
