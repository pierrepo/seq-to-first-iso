
"""Create a tsv file from a file of sequences.

Read a file composed of sequences of amino acids on each line and
return sequence, mass, formula, M0 and M1 in normal and C[12] conditions
as a tsv file.

Naming conventions for isotopes follow pyteomics's conventions.
"""

__authors__ = "Lilian Yang-crosson"
__copyright__ = ""
__credits__ = [""]
__license__ = ""
__version__ = "0.1"
__maintainer__ = ""
__email__ = ""

import argparse
from pathlib import Path

import pandas as pd
from pyteomics import mass


USAGE_ERROR = "Usage: python seq-to-first-iso.py filename [-o output] [-n aa]"
# Note: pyteomics also have H- and -OH that can be used for sequences
# which are not supported in this version.
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYUO")


def user_input():
    """
    Handle the user parameter from the command line.

    Parse the command line parameter and build the list of input files.
    """
    parser = argparse.ArgumentParser(
        description="Read a file of sequences and creates a tsv file")

    # Input file is required as a positional argument.
    parser.add_argument("input", type=Path, help="file to parse")

    # Optional arguments.
    parser.add_argument("-o", "--output", nargs="?", type=Path,
                        help="name of output file")
    parser.add_argument("-n", "--non_labelled_aa", nargs="+",
                        metavar="amino_a",
                        help="amino acids with default abundance")

    options = parser.parse_args()
    return options


def sequence_parser(file):
    """Return a list of peptide sequences parsed from a file.

    Take a Path() object as argument, return a list of uppercase peptides
    or exit if the list is empty.
    """
    # Check if file exists.
    if not file.is_file():
        exit("Error: file {} does not exist in current directory {}\n".format(
            file, file.cwd()) + USAGE_ERROR)

    # Obtain a list of sequences as string if they are amino acids.
    with open(file, "r") as filin:
        sequences = []
        ignored_lines = 0
        for sequence in filin:
            upper_sequence = sequence.upper()
            # Character not recognized as amino acid.
            if not set(upper_sequence) - AMINO_ACIDS:
                print(len(upper_sequence))
                sequences.append(upper_sequence)
            else:
                ignored_lines += 1

    # Check if the file format is correct.
    if not sequences:
        exit("Error: incorrect format, make sure that lines "
             "in {} are valid sequences of amino acids".format(file))
    if ignored_lines:
        print("{} lines ignored out of {}".format(ignored_lines,
                                                  ignored_lines+len(sequences)
                                                  ))
    return sequences


def compute_M0(f, a):
    """Return the monoisotopic abundance M0 of a sequence with its formula.

    f is the chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}
    a is the abundance of isotopes, as a dict in the format:
        {element_name[isotope_number]: relative abundance, ..}

    """
    M0 = a["C[12]"]**f["C"] * a["H[1]"]**f["H"] * a["N[14]"]**f["N"] \
        * a["O[16]"]**f["O"] * a["S[32]"]**f["S"]
    return M0


def compute_M1(f, a):
    """Compute abundance of second isotopologue M1 from its formula.

    f is the chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}
    a is the abundance of isotopes, as a dict in the format:
        {element_name[isotope_number]: relative abundance, ..}

    """
    M1 = (
          (f["C"] * a["C[12]"]**(f["C"]-1) * a["C[13]"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["H"] * a["C[12]"]**f["C"]
              * a["H[1]"]**(f["H"]-1) * a["H[2]"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["N"] * a["C[12]"]**f["C"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**(f["N"]-1) * a["N[15]"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["O"] * a["C[12]"]**f["C"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**(f["O"]-1) * a["O[17]"]
              * a["S[32]"]**f["S"])

          + (f["S"] * a["C[12]"]**f["C"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**(f["S"]-1) * a["S[33]"])
          )
    return M1


if __name__ == "__main__":
    # Default isotopic abundances from MIDAs website:
    # https://www.ncbi.nlm.nih.gov/CBBresearch/Yu/midas/index.html .
    isotopic_abundance = {"H[1]": 0.999885, "H[2]": 0.000115,
                          "C[12]": 0.9893,  "C[13]": 0.0107,
                          "N[14]": 0.99632, "N[15]": 0.00368,
                          "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                          "S[32]": 0.9493,  "S[33]": 0.0076,  "S[34]": 0.0429}

    options = user_input()
    input_file = options.input
    auxotrophic_aa = options.non_labelled_aa
    print(auxotrophic_aa)
    sequences = sequence_parser(input_file)

    # Dataframe of sequences.
    df_peptides = pd.DataFrame({"sequence": sequences})

    # Add mass and formulas of sequences
    df_peptides["mass"] = df_peptides["sequence"].map(mass.calculate_mass)
    df_peptides["formula"] = df_peptides["sequence"].map(mass.Composition)

    # Add M0 and M1 in normal conditions.
    df_peptides["M0_NC"] = df_peptides["formula"].apply(compute_M0,
                                                        a=isotopic_abundance)
    df_peptides["M1_NC"] = df_peptides["formula"].apply(compute_M1,
                                                        a=isotopic_abundance)

    # Add M0 and M1 in 99.99 % of C[12].
    prop = 0.9999
    isotopic_abundance["C[12]"] = prop
    isotopic_abundance["C[13]"] = 1-prop
    df_peptides["M0_12C"] = df_peptides["formula"].apply(compute_M0,
                                                         a=isotopic_abundance)
    df_peptides["M1_12C"] = df_peptides["formula"].apply(compute_M1,
                                                         a=isotopic_abundance)

    # Choose output filename.
    if not options.output:
        output_file = input_file.stem + ".tsv"
    else:
        output_file = options.output

    # Import dataframe to tsv file.
    df_peptides.to_csv(output_file, sep="\t")
