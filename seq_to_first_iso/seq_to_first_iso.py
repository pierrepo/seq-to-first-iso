
"""Create a tsv file from a file of sequences.

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

import argparse
import logging
from pathlib import Path
import sys

import pandas as pd
from pyteomics import mass


USAGE_ERROR = "Usage: python seq-to-first-iso.py filename " \
             + "[-o output] [-n aa]"
# Note: pyteomics also have U, O, H- and -OH that can be used for sequences
# which are not supported in this version.
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

# Set custom logger.
log = logging.getLogger(__name__)
log_handler = logging.StreamHandler()
log_formatter = logging.Formatter("[%(asctime)s] %(levelname)-8s: %(message)s",
                                  "%Y-%m-%d, %H:%M:%S")
log_handler.setFormatter(log_formatter)
log.addHandler(log_handler)
log.setLevel(logging.INFO)


# Default isotopic abundances from MIDAs website:
# https://www.ncbi.nlm.nih.gov/CBBresearch/Yu/midas/index.html .
# X is C with default abundance.
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


def user_input(args):
    """Parse and handle the submitted command line."""
    parser = argparse.ArgumentParser(
        description="Read a file of sequences and creates a tsv file")

    # Input file is required as a positional argument.
    parser.add_argument("input", type=Path, help="file to parse")

    # Optional arguments.
    parser.add_argument("-o", "--output", type=str,
                        help="name of output file")
    parser.add_argument("-n", "--non-labelled-aa",
                        metavar="amino_a",
                        help="amino acids with default abundance")

    options = parser.parse_args(args)

    # Check if file exists.
    if not options.input.is_file():
        log.error(f"file {options.input} does not exist in "
                  f"current directory '{options.input.cwd()}'\n{USAGE_ERROR}"
                  )
        exit()

    # Check if amino acids are correct.  If not, tell which one.
    if not options.non_labelled_aa:
        # Change to empty list to avoid Nonetype errors.
        options.non_labelled_aa = []
    else:
        options.non_labelled_aa = options.non_labelled_aa.split(",")
        # Convert amino acids to uppercase for compatibility.
        options.non_labelled_aa = [char.upper()
                                   for char in options.non_labelled_aa]
        unrecognized_aa = []

        for arg in options.non_labelled_aa:
            if arg not in AMINO_ACIDS:
                unrecognized_aa.append(arg)

        if unrecognized_aa:
            log.warning(f"{unrecognized_aa} not recognized as amino acid")

    return options


def sequence_parser(file):
    """Return a tuple (sequences:list, ignored_lines:int) parsed from a file.

    Take a Path() object as argument, return a list of uppercase peptides
    and the number of ignored lines.
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

    return sequences, ignored_lines


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


def seq_to_midas(sequence_l, sequence_nl):
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


def seq_to_tsv(sequences, output_file, unlabelled_aa):
    """Create a tsv from sequences and return its name.

    Take a list of amino acid sequences, a string for the output filename
    and a list of unlabelled amino acids.
    """
    # Dataframe of sequences.
    df_peptides = pd.DataFrame({"sequence": sequences})

    # Separate sequences.
    df_peptides["labelled"], df_peptides["unlabelled"] = zip(
            *df_peptides["sequence"].apply(separate_labelled,
                                           unlabelled_aa=unlabelled_aa))

    # Add mass and formulas of sequences
    log.info("Computing mass")
    df_peptides["mass"] = df_peptides["sequence"].map(mass.calculate_mass)

    log.info("Computing formula")
    # Formula as a string (instead of mass.Composition).
    df_peptides["f"] = df_peptides["sequence"].apply(mass.Composition)
    df_peptides["formula"] = df_peptides["f"].apply(formula_to_str)
    # Composition, with unlabelled C as element X.
    df_peptides["f_X"] = df_peptides.apply(lambda x:
                                           seq_to_midas(x["labelled"],
                                                        x["unlabelled"]),
                                           axis=1)

    # Add M0 and M1 in normal conditions.
    log.info("Computing M0 and M1")
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

    wanted_columns = ["sequence", "mass", "formula", "formula_X",
                      "M0_NC", "M1_NC", "M0_12C", "M1_12C"]
    # Import dataframe to tsv file.
    #df_peptides[wanted_columns].to_csv(output_file, sep="\t", index=False)

    return df_peptides[wanted_columns]


def cli(args=None):
    """Entry point for seq_to_first_iso's CLI."""
    if not args:
        args = sys.argv[1:]

    options = user_input(args)
    input_file = options.input
    unlabelled_aa = options.non_labelled_aa

    if unlabelled_aa:
        log.info(f"Amino acid with default abundance: {unlabelled_aa}")

    log.info("Parsing file")
    sequences, ignored_lines = sequence_parser(input_file)

    if not sequences:
        log.error(f"incorrect format, make sure that lines "
                  f"in {str(input_file)} are valid sequences of amino acids")
        exit()
    if ignored_lines:
        log.warning(f"{ignored_lines} lines ignored out of "
                    f"{ignored_lines+len(sequences)}")

    # Choose output filename.
    if not options.output:
        output_file = input_file.stem + ".tsv"
    else:
        output_file = options.output + ".tsv"

    seq_to_tsv(sequences, output_file, unlabelled_aa)


if __name__ == "__main__":
    cli()
