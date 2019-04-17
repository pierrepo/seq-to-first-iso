
"""Seq to first iso: compute first isotopologues intensities from sequences.

Create a tsv file from a file of sequences.

Read a file composed of sequences of amino acids on each line and return :
    sequence,
    mass,
    normal formula and formula with unlabelled amino acids,
    M0 and M1 in normal and C[12] conditions
as a tsv file.

Unlabelled amino acids's carbons keep default isotopic abundance,
and are represented as X in formulas.
Naming conventions for isotopes follow pyteomics's conventions.
"""

__authors__ = "Lilian Yang-crosson, Pierre Poulain"
__license__ = "BSD 3-Clause License"
__version__ = "0.2.0"
__maintainer__ = "Pierre Poulain"
__email__ = "pierre.poulain@cupnet.net"

from .seq_to_first_iso import (AMINO_ACIDS,
                               C12_abundance,
                               isotopic_abundance,
                               sequence_parser,
                               separate_labelled,
                               compute_M0_nl,
                               compute_M1_nl,
                               formula_to_str,
                               seq_to_midas,
                               seq_to_tsv,
                               )
__all__ = ["AMINO_ACIDS",
           "C12_abundance",
           "isotopic_abundance",
           "sequence_parser",
           "separate_labelled",
           "compute_M0_nl",
           "compute_M1_nl",
           "formula_to_str",
           "seq_to_midas",
           "seq_to_tsv",
           ]
