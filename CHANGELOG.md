# Changelog

**Dev**

**0.3.0** (2019-04-18)

*Changed*
- Add support for files with annotations before the sequences
  - sequence_parser() now returns (annotations, sequences, ignored_lines)
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
