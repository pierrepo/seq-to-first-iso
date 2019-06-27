
"""Compute first two isotopologue intensities from sequences.

The program computes M0 and M1 and differentiate labelled
(with a 99.99 % C[12] enrichment) and unlabelled amino acids.

Read a file composed of amino acid sequences on each line and return :
    sequence, mass, formula, formula_X, M0_NC, M1_NC, M0_12C and M1_12C
in a tsv file.
- formula_X is the chemical formula with carbon of unlabelled
amino acids marked as X.
- NC are Normal Conditions, 12C are C[12] enrichment Conditions.


Example
-------
Running the script after installation
    $ seq-to-first-iso sequences.txt
will provide file 'sequences_stfi.tsv'


Notes
-----
Carbon of unlabelled amino acids keep default isotopic abundance,
and are represented as X in formulas.
Naming conventions for isotopes follow pyteomics's conventions.

"""

__authors__ = "Lilian Yang-crosson, Pierre Poulain"
__license__ = "BSD 3-Clause License"
__version__ = "0.4.3"
__maintainer__ = "Pierre Poulain"
__email__ = "pierre.poulain@cupnet.net"

from .seq_to_first_iso import (AMINO_ACIDS,
                               C12_abundance,
                               isotopic_abundance,
                               UNIMOD_MODS,
                               sequence_parser,
                               separate_labelled,
                               compute_M0_nl,
                               compute_M1_nl,
                               formula_to_str,
                               seq_to_xcomp,
                               get_mods_composition,
                               seq_to_df,
                               )
__all__ = ["AMINO_ACIDS",
           "C12_abundance",
           "isotopic_abundance",
           "UNIMOD_MODS",
           "sequence_parser",
           "separate_labelled",
           "compute_M0_nl",
           "compute_M1_nl",
           "formula_to_str",
           "seq_to_xcomp",
           "get_mods_composition",
           "seq_to_df",
           ]
