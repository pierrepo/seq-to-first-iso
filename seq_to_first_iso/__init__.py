"""Compute intensities of the first two isotopologue.

Use peptide sequences and charges.

The program computes M0 and M1 and differentiate labelled
(with a 99.99 % C[12] enrichment) and unlabelled amino acids.

Read a .tsv file composed of amino acid sequences on each line and return:
sequence, mass, formula, formula_X, M0_NC, M1_NC, M0_12C and M1_12C
in a .tsv file.

Formula_X is the chemical formula with carbon of unlabelled
amino acids marked as X.

NC means Normal Condition, 12C means C[12] enrichment condition.


Example
-------
Running the script after installation

    $ seq-to-first-iso sequences.tsv sequence_column_name charge_column_name

will provide file 'sequences_stfi.tsv'


Notes
-----
Carbon of unlabelled amino acids keep default isotopic abundance,
and are represented as X in formulas.
Naming conventions for isotopes follow pyteomics's conventions.

"""

__authors__ = "Lilian Yang-crosson, Pierre Poulain"
__license__ = "BSD 3-Clause License"
__version__ = "1.1.0"
__maintainer__ = "Pierre Poulain"
__email__ = "pierre.poulain@cupnet.net"

from .seq_to_first_iso import (AMINO_ACIDS,
                               XTANDEM_MOD_PATTERN,
                               UNIMOD_MODS,
                               NATURAL_ABUNDANCE,
                               C12_ABUNDANCE,
                               parse_input_file,
                               filter_input_dataframe,
                               check_amino_acids,
                               separate_labelled,
                               compute_M0_nl,
                               compute_M1_nl,
                               formula_to_str,
                               convert_atom_C_to_X,
                               get_charge_composition,
                               get_mods_composition,
                               compute_intensities,
                               )

__all__ = ["AMINO_ACIDS",
           "XTANDEM_MOD_PATTERN",
           "UNIMOD_MODS",
           "NATURAL_ABUNDANCE",
           "C12_ABUNDANCE",
           "parse_input_file",
           "filter_input_dataframe",
           "check_amino_acids",
           "separate_labelled",
           "compute_M0_nl",
           "compute_M1_nl",
           "formula_to_str",
           "convert_atom_C_to_X",
           "get_charge_composition",
           "get_mods_composition",
           "compute_intensities",
           ]
