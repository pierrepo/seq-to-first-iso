[![PyPI version](https://badge.fury.io/py/seq-to-first-iso.svg)](https://badge.fury.io/py/seq-to-first-iso)
![Conda](https://img.shields.io/conda/v/bioconda/seq-to-first-iso.svg?color=bright-green)
[![Build Status](https://travis-ci.org/pierrepo/seq-to-first-iso.svg?branch=master)](https://travis-ci.org/pierrepo/seq-to-first-iso)
[![Documentation Status](https://readthedocs.org/projects/seq-to-first-iso/badge/?version=latest)](https://seq-to-first-iso.readthedocs.io/en/latest/?badge=latest)

# Seq-to-first-iso

> Compute first two isotopologues intensity from peptide sequence

The program differentiate labelled and unlabelled amino acids
for a 99.99 % 12C enrichment.  
The documentation can be found [here](https://seq-to-first-iso.readthedocs.io/en/latest/).  
Try the demo with Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/pierrepo/seq-to-first-iso/master)

## Installation

### With pip

```
$ pip install seq-to-first-iso
```

### With conda

```
$ conda install seq-to-first-iso -c bioconda
```

### Developer mode


Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

Clone repo:
```
$ git clone https://github.com/pierrepo/seq-to-first-iso
$ cd seq-to-first-iso
```

Create conda environment:
```
$ conda env create -f environment.yml
```

Remark: for a fully reproducible environment, you could also use:
```
$ conda env create -f environment.lock.yml
```


Activate conda environment:
```
$ conda activate seq-to-first-iso
```

Install local package:
```
$ pip install -e .
```

## Usage

The script takes a file with one sequence of amino acids per line and returns a tsv of the file with columns:

|sequence|mass|formula|formula_X| M0_NC | M1_NC | M0_12C | M1_12C |
|--------|----|-------|---------|-------|-------|--------|--------|

Once installed, the script can be called with:

```shell
$ seq-to-first-iso filename [-o output_name] [-n amino_acids...]
```
Optional arguments are in square brackets  
This will create _filename_stfi.tsv_ if filename is a correct file

**0.3.0 :** The input file can have annotations separated by a tabulation before the sequences  
**0.4.0 :** Support for [X!Tandem](https://www.thegpm.org/tandem/) Post-Translational Modifications added

### Options

- `-h, --help`:  
Provide a help page

- `-v, --version`:  
Provide the version

- `-o, --output`:  
Change the name of the output file

- `-n, --non-labelled-aa`:  
Take 1 or more amino acid separated by a comma


### Examples

- You can provide a list of amino acids which will keep default isotopic abundance:

Supposing [peptides.txt](https://github.com/pierrepo/seq-to-first-iso/blob/master/notebooks/peptides.txt) :

```
YAQEISR
VGFPVLSVKEHK
LAMVIIKEFVDDLK
```

The command
```shell
$ seq-to-first-iso peptides.txt -n V,W
```
will create [peptides_stfi.tsv](https://github.com/pierrepo/seq-to-first-iso/blob/master/notebooks/peptides_stfi.tsv) :

|sequence| mass| formula|formula_X| M0_NC| M1_NC| M0_12C| M1_12C|
|--------|-----|--------|---------|------|------|-------|-------|
YAQEISR| 865.42938099921| C37H59O13N11| C37H59O13N11| 0.6206414140575179|	0.280870823368276| 0.9206561231798033| 0.05161907174495234|
VGFPVLSVKEHK| 1338.7659712609| C63H102O16N16| C48H102O16N16X15|  0.4550358985377136| 0.34506032928190855| 0.7589558393662944| 0.18515489894512063|
LAMVIIKEFVDDLK| 1632.91606619252| C76H128O21N16S1| C66H128O21N16S1X10| 0.36994021481230627| 0.3373188347614264| 0.7475090558698947| 0.15292723586285323|

Where, in 12C enrichment conditions, the isotopologue intensity M0_12C and M1_12C are computed with unlabelled Valine and Tryptophan (V and W have default isotopic abundance)


- You can change the name of the output file:

```shell
$ seq-to-first-iso peptides.txt -o sequence
```
will create a file named *sequence.tsv*


## Credits

- **Binder**
  - Jupyter et al., "Binder 2.0 - Reproducible, Interactive, Sharable
  Environments for Science at Scale." Proceedings of the 17th Python
  in Science Conference. 2018. 10.25080/Majora-4af1f417-011

- **Bioconda**:
  - Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. “Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences”. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

- **MIDAs**:
  - Alves G, Ogurtsov AY, Yu YK (2014) Molecular Isotopic Distribution Analysis (MIDAs) with adjustable mass accuracy. J Am Soc Mass Spectrom, 25: 57-70. DOI: 10.1007/s13361-013-0733-7

- **Pyteomics**:
  - Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

  - Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
