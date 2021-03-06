# Changelog

**Dev**

**1.1.0** (2019-12-19)
- Add export_to_knime() function

**1.0.0** (2019-12-05)
- Simplify interface
- Take .tsv as input file
- Take peptide sequence and charge as input data (from input file)
- Check input file and input dataframe
- Comply with PEP 8 and PEP 257
- Update API and CLI demo notebooks
  
**0.5.1** (2019-07-10)

*Fixed*
- seq_to_xcomp() can now correctly take a pyteomics.mass.Composition object as the second parameter

**0.5.0** (2019-07-03)

*Added*
- Flag for version number in CLI
- Jupyter notebook with Binder environment for demonstrations
- Documentation on Read the Docs
- Conda availability

*Changed*
- **Breaking change** : changed seq_to_midas() to seq_to_xcomp()
- **Breaking change** : changed seq_to_tsv() to seq_to_df()

**0.4.3** (2019-06-26)

*Fixed*
- Fix requirements not being installed with `pip install`

**0.4.2** (2019-06-25)

*Fixed*
- Fix *setup.cfg*'s installation requirements

**0.4.1** (2019-06-24)

- Extend numpydoc style to all functions in *seq_to_first_iso.py*

**0.4.0** (2019-06-21)

*Changed*
- Add support for Xtandem Parsing
  - **Breaking change:** sequence_parser() now returns a dict with "annotations", "raw_sequences", "sequences", "modifications" and "ignored_lines"
  - Add get_mods_composition() that returns a composition from a list of Unimod PTMs
  - Remove the appended "\_stfi" if -o flag is provided

**0.3.0** (2019-04-18)

*Changed*
- Add support for files with annotations before the sequences
  - **Breaking change:** sequence_parser() now returns (annotations, sequences, ignored_lines)
  - seq_to_tsv() now accepts (sequences, unlabelled_aa, annotations=None)

*Fixed*
- Output files now have "\_stfi" appended to differentiate from .tsv input files with the same name

**0.2.1** (2019-04-17)
- Format CHANGELOG

**0.2.0** (2019-04-17)
- Add bumpversion

*Changed*
- seq_to_tsv() no longer writes a file, instead it returns a dataframe

**0.1.0** (2019-04-08)
- First release
