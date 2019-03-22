
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


USAGE_ERROR = "Usage: python seq-to-first-iso.py filename " \
             + "[-o output] [-n aa]"
# Note: pyteomics also have H- and -OH that can be used for sequences
# which are not supported in this version.  They are implicitly added.
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

    # Check if file exists.
    if not options.input.is_file():
        exit("Error: file {} does not exist in current directory {}\n".format(
            options.input, options.input.cwd()) + USAGE_ERROR)

    # Check if amino acids are correct.  If not, tell which one.
    unlabelled_aa = options.non_labelled_aa
    if not unlabelled_aa:
        # Change to empty list to avoid Nonetype errors.
        options.non_labelled_aa = []
    elif unlabelled_aa:
        # Convert amino acids to uppercase for compatibility.
        options.non_labelled_aa = [char.upper() for char in unlabelled_aa]
        unrecognized_aa = []

        for arg in options.non_labelled_aa:
            if arg not in AMINO_ACIDS:
                unrecognized_aa.append(arg)

        if unrecognized_aa:
            print("Warning: {} not recognized "
                  "as amino acid".format(unrecognized_aa))

    return options


def sequence_parser(file):
    """Return a list of peptide sequences parsed from a file.

    Take a Path() object as argument, return a list of uppercase peptides
    or exit if the list is empty.
    """
    # Obtain a list of sequences as string if they are amino acids.
    with open(file, "r") as filin:
        sequences = []
        ignored_lines = 0
        for sequence in filin:
            upper_sequence = sequence.upper().strip()
            # Character not recognized as amino acid.
            if not (set(upper_sequence) - AMINO_ACIDS) and upper_sequence:
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


def separate_labelled(sequence, unlabelled_aa):
    """Get the sequence of unlabelled amino acids from a sequence.

    sequence is a string of amino acids
    unlabelled_aa is a container (list, string...) of unlabelled amino acids

    Return the sequences as a tuple of string with:
        the sequence without the unlabelled amino acids
        the unlabelled amino acids in the sequence

    Note: we can also use comprehension lists (might go faster)
    labelled_seq = "".join([char in sequence if char not in unlabelled_aa])
    unlabelled_seq = "".join([char in sequence if char in unlabelled_aa])
    """
    labelled_seq = []
    unlabelled_seq = []
    for char in sequence:
        if char in unlabelled_aa:
            unlabelled_seq.append(char)
        else:
            labelled_seq.append(char)
    return "".join(labelled_seq), "".join(unlabelled_seq)


def compute_M0_nl(f, a):
    """Return the monoisotopic abundance M0 of a formula with mixed labels.

    f is the chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}
    a is the abundance of isotopes, as a dict in the format:
        {element_name[isotope_number]: relative abundance, ..}

    X represents C with default isotopic abundance.
    """
    M0 = a["C[12]"]**f["C"] * a["X[12]"]**f["X"] * a["H[1]"]**f["H"] \
        * a["N[14]"]**f["N"] * a["O[16]"]**f["O"] * a["S[32]"]**f["S"]
    return M0


def compute_M1_nl(f, a):
    """Compute abundance of second isotopologue M1 from its formula.

    f is the chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}
    a is the abundance of isotopes, as a dict in the format:
        {element_name[isotope_number]: relative abundance, ..}

    X represents C with default isotopic abundance.
    """
    M1 = (
          (f["C"] * a["C[12]"]**(f["C"]-1) * a["C[13]"]
              * a["X[12]"]**f["X"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["X"] * a["C[12]"]**f["C"]
              * a["X[12]"]**(f["X"]-1) * a["X[13]"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["H"] * a["C[12]"]**f["C"]
              * a["X[12]"]**f["X"]
              * a["H[1]"]**(f["H"]-1) * a["H[2]"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["N"] * a["C[12]"]**f["C"]
              * a["X[12]"]**f["X"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**(f["N"]-1) * a["N[15]"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**f["S"])

          + (f["O"] * a["C[12]"]**f["C"]
              * a["X[12]"]**f["X"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**(f["O"]-1) * a["O[17]"]
              * a["S[32]"]**f["S"])

          + (f["S"] * a["C[12]"]**f["C"]
              * a["X[12]"]**f["X"]
              * a["H[1]"]**f["H"]
              * a["N[14]"]**f["N"]
              * a["O[16]"]**f["O"]
              * a["S[32]"]**(f["S"]-1) * a["S[33]"])
          )

    return M1


def formula_to_str(composition):
    """Return formula from Composition as a string."""
    formula_str = ""
    for element in "CHONSX":
        if element in composition:
            formula_str += f"{element}{composition[element]}"
    return formula_str


def seq_to_midas(sequence_l, sequence_nl):  # ADDED
    """Take 2 amino acid sequences and return the formula for MIDAs.

    Note: the function assumes the second sequence has no terminii.
    """
    formula_l = mass.Composition(sequence_l)
    formula_nl = mass.Composition(parsed_sequence=sequence_nl)
    try:
        formula_nl["X"] = formula_nl.pop("C")
    except KeyError:
        pass
    return formula_l+formula_nl


def formula_to_midas(formula_l, formula_nl):
    """Take 2 chemical formulas and return the formula for MIDAs.

    Note: can adapt to take sequence to composition.
    """
    f2 = mass.Composition(formula_nl)  # to not change formula2
    try:
        f2["X"] = f2.pop("C")  # Based on the assumption there are C.
    except KeyError:
        pass
    return formula_l+f2


if __name__ == "__main__":
    # Default isotopic abundances from MIDAs website:
    # https://www.ncbi.nlm.nih.gov/CBBresearch/Yu/midas/index.html .
    # X is C with default abundance
    isotopic_abundance = {"H[1]": 0.999885, "H[2]": 0.000115,
                          "C[12]": 0.9893,  "C[13]": 0.0107,
                          "X[12]": 0.9893,  "X[13]": 0.0107,
                          "N[14]": 0.99632, "N[15]": 0.00368,
                          "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                          "S[32]": 0.9493,  "S[33]": 0.0076,  "S[34]": 0.0429}

    C12_abundance = dict(isotopic_abundance)
    prop = 0.9999
    C12_abundance["C[12]"] = prop
    C12_abundance["C[13]"] = 1-prop

    options = user_input()
    input_file = options.input
    unlabelled_aa = options.non_labelled_aa

    if unlabelled_aa:
        print("Amino acid with default abundance: {}".format(unlabelled_aa))

    print("Parsing file")
    sequences = sequence_parser(input_file)

    # Dataframe of sequences.
    df_peptides = pd.DataFrame({"sequence": sequences})

    # Separate sequences.
    df_peptides["labelled"], df_peptides["unlabelled"] = zip(
            *df_peptides["sequence"].apply(separate_labelled,
                                           unlabelled_aa=unlabelled_aa))

    # Add mass and formulas of sequences
    print("Computing mass")
    df_peptides["mass"] = df_peptides["sequence"].map(mass.calculate_mass)

    print("Computing formula")
    # Formula as a string (instead of mass.Composition).
    df_peptides["f"] = df_peptides["sequence"].apply(mass.Composition)
    df_peptides["formula"] = df_peptides["f"].apply(formula_to_str)
    # Composition, with unlabelled C as element X.
    df_peptides["f_X"] = df_peptides.apply(lambda x:
                                           seq_to_midas(x["labelled"],
                                                        x["unlabelled"]),
                                           axis=1)

    # Add M0 and M1 in normal conditions.
    print("Computing M0 and M1")
    # Can use compute_M0_nl with isotopic abundance twice
    df_peptides["M0_NC"] = df_peptides["f_X"].apply(compute_M0_nl,
                                                    a=isotopic_abundance)
    df_peptides["M1_NC"] = df_peptides["f_X"].apply(compute_M1_nl,
                                                    a=isotopic_abundance)

    df_peptides["M0_12C"] = df_peptides["f_X"].apply(compute_M0_nl,
                                                     a=C12_abundance)
    df_peptides["M1_12C"] = df_peptides["f_X"].apply(compute_M1_nl,
                                                     a=C12_abundance)

    # For verification with MIDAs, might be removed.
    df_peptides["formula_X"] = df_peptides["f_X"].apply(formula_to_str)

    # Choose output filename.
    if not options.output:
        output_file = input_file.stem + ".tsv"
    else:
        output_file = options.output

    wanted_columns = ["sequence", "mass", "formula", "formula_X",
                      "M0_NC", "M1_NC", "M0_12C", "M1_12C"]
    # Import dataframe to tsv file.
    df_peptides[wanted_columns].to_csv(output_file, sep="\t", index=False)
