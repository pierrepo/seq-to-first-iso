"""Main module of seq_to_first_iso.

Provide function to compute M0 and M1 with labelled and unlabelled amino acids
for the case of a 99.99 % C[12] enrichment.
The Command-Line Interface is defined here.


Example
-------
Running the script directly
    $ python seq_to_first_iso sequences.txt
will provide file 'sequences_stfi.tsv'


Notes
-----
Carbon of unlabelled amino acids keep default isotopic abundance,
and are represented as X in formulas.
Naming conventions for isotopes follow pyteomics's conventions.


Attributes
----------
USAGE_ERROR : str
    Default message for errors.
AMINO_ACIDS : set
    Set of supported 1-letter amino acids.
XTANDEM_MOD_PATTERN : re.Pattern
    Regular expression capturing XTandem Post Translational Modifications.
UNIMOD_MODS : pyteomics.mass.Unimod
    Dictionary with Unimods entries.
USED_ELEMS : str
    String of elements used/recognized by the program.
isotopic_abundance : dict
    Dictionary of isotopic abundances with values taken from MIDAs.
C12_abundance : dict
    Dictionary of isotopic abundances with C[12] abundance at 0.9999.
log : logging.Logger
    Logger outputting in text terminals.

"""

import argparse
import logging
from pathlib import Path
import re
import sys

import pandas as pd
from pyteomics import mass

from seq_to_first_iso import __version__


USAGE_ERROR = "Usage: python seq-to-first-iso.py filename.tsv sequence_column_name charge_column_name " \
             + "[-o output] [-u aa]"
# Note: pyteomics also have U, O, H- and -OH that can be used for sequences
# which are not supported in this version.
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

XTANDEM_MOD_PATTERN = re.compile(r"""
                                 \.?       # 0 or 1 dot
                                 \(        # Opening parenthesis
                                   (       # Begin capture
                                    (?:        # Not capture the following
                                     [^\(\)] | # Either not parentheses or
                                     \(-?\d+\) # parentheses containing an int
                                    )+         # multiple times
                                   )       # End capture
                                 \)        # Closing parenthesis
                                 """, re.VERBOSE)

UNIMOD_MODS = mass.Unimod()
# This variable is obsoleted if an natural element shall be named X.
USED_ELEMS = "CHONPSX"

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
normal_abundance = {"H[1]": 0.999885, "H[2]": 0.000115,
                    "C[12]": 0.9893,  "C[13]": 0.0107,
                    "X[12]": 0.9893,  "X[13]": 0.0107,
                    "N[14]": 0.99632, "N[15]": 0.00368,
                    "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                    "S[32]": 0.9493,  "S[33]": 0.0076,  "S[34]": 0.0429}

C12_abundance = dict(normal_abundance)
prop = 0.9999
C12_abundance["C[12]"] = prop
C12_abundance["C[13]"] = 1-prop


def user_input(args):
    """Parse and handle the submitted command line.

    Parameters
    ----------
    args : list of str
        List of arguments received from the CLI.

    Returns
    -------
    argparse.Namespace
        Object containing the arguments parsed from the CLI.

    Raises
    ------
    SystemExit
        If the file provided is not found.

    """
    parser = argparse.ArgumentParser(
        description="Read a tsv file with sequences and charges and compute intensity of first isotopologues")

    # Input file is required as a positional argument.
    parser.add_argument("input_file_name", type=Path, help="file to parse in .tsv format")
    parser.add_argument("sequence_col_name", type=str, help="column name with sequences")
    parser.add_argument("charge_col_name", type=str, help="column name with charges")

    # Optional arguments.
    parser.add_argument("-o", "--output", type=str,
                        help="name of output file")
    parser.add_argument("-u", "--unlabelled-aa",
                        metavar="amino_a",
                        help="amino acids with default abundance")

    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(__version__))

    options = parser.parse_args(args)

    # Check if file exists.
    if not options.input_file_name.is_file():
        log.error(f"file {options.input_file_name} does not exist in "
                  f"current directory '{options.input_file_name.cwd()}'\n{USAGE_ERROR}"
                  )
        exit()

    # Check if amino acids are correct.  If not, tell which one.
    if not options.unlabelled_aa:
        # Change to empty list to avoid Nonetype errors.
        options.unlabelled_aa = []
    else:
        options.unlabelled_aa = options.unlabelled_aa.split(",")
        # Convert amino acids to uppercase for compatibility.
        options.unlabelled_aa = [char.upper()
                                   for char in options.unlabelled_aa]
        unrecognized_aa = []

        for arg in options.unlabelled_aa:
            if arg not in AMINO_ACIDS:
                unrecognized_aa.append(arg)

        if unrecognized_aa:
            log.warning(f"{unrecognized_aa} not recognized as amino acid")

    return options


def parse_input_file(filename, sequence_col_name, charge_col_name, sep="\t"):
    r"""Parse input file with peptide sequences and charges.

    Parameters
    ----------
    filename : str
        Filename, the file can either just have sequences for each line or
        can have have annotations and sequences with a separator in-between.
    sequence_col_name : str
        Name of column with peptide sequences
    charge_col_name : str
        Name of column with peptide charges
    sep : str, optional
        Separator for files with annotations (default is ``\t``).

    Returns
    -------
    pandas dataframe
        | With columns :
        |     - "sequences": peptide sequences.
        |     - "charges": peptide charges.

    """
    df = pd.read_csv(filename, sep='\t')
    line_count, row_count = df.shape
    log.info(f"Read {filename} with {line_count} lines and {row_count} columns")
    if sequence_col_name not in df.columns:
        log.error(f"Column '{sequence_col_name}' not found in data.")
    if charge_col_name not in df.columns:
        log.error(f"Column '{charge_col_name}' not found in data.")
    # Keep only sequences and charges column from original dataframe.
    df = df[[sequence_col_name, charge_col_name]]
    # Rename sequences and charges column to internal naming scheme.
    df.rename(columns={sequence_col_name: "sequence", charge_col_name: "charge"})
    return df


def check_amino_acids(seq):
    r""" Check elements of a sequence are known amino acids.

    Parameters
    ----------
    seq : str
        Peptide sequence.
    
    Returns
    -------
    Tuple of two str
        | (sequence, "") if the sequence is composed 
        |                of allowed amino acids
        | ("", "Unrecognized amino acids") if the sequence is composed 
        |                                  of unallowed amino acids.
    """
    if not(set(seq) - AMINO_ACIDS) and seq:
        return seq, ""
    else:
        return "", "Unrecognized amino acids"


def compute_M0(f, a):
    """Return the monoisotopic abundance M0 of a sequence with its formula.

    Parameters
    ----------
    f : pyteomics.mass.Composition
        Chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}.
    a : dict
        Dictionary of abundances of isotopes, in the format:
        {element_name[isotope_number]: relative abundance, ..}.

    Returns
    -------
    float
        Value of M0.

    Notes
    -----
    Unused. Use compute_M0_nl instead.

    """
    M0 = a["C[12]"]**f["C"] * a["H[1]"]**f["H"] * a["N[14]"]**f["N"] \
        * a["O[16]"]**f["O"] * a["S[32]"]**f["S"]
    return M0


def compute_M1(f, a):
    """Compute abundance of second isotopologue M1 from its formula.

    Parameters
    ----------
    f : pyteomics.mass.Composition
        Chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}.
    a : dict
        Dictionary of abundances of isotopes, in the format:
        {element_name[isotope_number]: relative abundance, ..}.

    Returns
    -------
    float
        Value of M1.

    Notes
    -----
    Unused. Use compute_M1_nl instead.

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

    Parameters
    ----------
    sequence : str
        String of amino acids.
    unlabelled_aa : container object
        Container (list, string...) of unlabelled amino acids.

    Returns
    -------
    tuple(str, str)
        | The sequences as a tuple of string with:
        |    - the sequence without the unlabelled amino acids
        |    - the unlabelled amino acids in the sequence

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

    Parameters
    ----------
    f : pyteomics.mass.Composition
        Chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}.
    a : dict
        Dictionary of abundances of isotopes, in the format:
        {element_name[isotope_number]: relative abundance, ..}.

    Returns
    -------
    float
        Value of M0.

    Notes
    -----
    X represents C with default isotopic abundance.

    """
    M0 = a["C[12]"]**f["C"] * a["X[12]"]**f["X"] * a["H[1]"]**f["H"] \
        * a["N[14]"]**f["N"] * a["O[16]"]**f["O"] * a["S[32]"]**f["S"]
    return M0


def compute_M1_nl(f, a):
    """Compute abundance of second isotopologue M1 from its formula.

    Parameters
    ----------
    f : pyteomics.mass.Composition
        Chemical formula, as a dict of counts for each element:
        {element_name: count_of_element_in_sequence, ...}.
    a : dict
        Dictionary of abundances of isotopes, in the format:
        {element_name[isotope_number]: relative abundance, ..}.

    Returns
    -------
    float
        Value of M1.

    Notes
    -----
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
    """Return formula from Composition as a string.

    Parameters
    ----------
    composition : pyteomics.mass.Composition
        Chemical formula.

    Returns
    -------
    str
        Human-readable string of the formula.

    Warnings
    --------
    If the composition has elements not in USED_ELEMS, they will not
    be added to the output.

    """
    formula_str = ""
    for element in USED_ELEMS:
        if element in composition:
            formula_str += f"{element}{composition[element]}"
    return formula_str


def convert_atom_C_to_X(sequence):
    """Replace carbon atom by element X atom in a composition.

    Parameters
    ----------
    sequence : str or pyteomics.mass.Composition
        Sequence or composition.

    Returns
    -------
    pyteomics.mass.Composition
        Composition with carbon atoms replaced by element X atoms.

    """
    # Force input to be a pyteomics.mass.Composition object.
    formula = mass.Composition(sequence)
    # Replace C atoms by X atoms.
    formula["X"] = formula.pop("C", 0)
    return formula


def get_charge_composition(charge):
    """Return the composition of a given charge (only H+).

    Parameters
    ----------
    charge : int
        Peptide charge.
    
    Returns
    -------
    pyteomics.mass.Composition
        Composition of the change (H+).
    """
    charge_composition = mass.Composition()
    charge_composition["H"] = charge
    return charge_composition


def get_mods_composition(modifications):
    """Return the composition of a list of modifications.

    Parameters
    ----------
    modifications : list of str
        List of modifications string (corresponding to Unimod titles).

    Returns
    -------
    pyteomics.mass.Composition
        The total composition change.

    """
    # ???: Have the mass.Unimod() dict as parameter ?
    total_mod_composition = mass.Composition()
    for mod in modifications:
        try:
            mod_composition = UNIMOD_MODS.by_title(mod)["composition"]
            total_mod_composition += mod_composition
            # Using set comparison here won't work with elements as isotopes.
            for elem in mod_composition:
                if elem not in USED_ELEMS:
                    log.warning(f"{elem} in ({mod}) is not supported "
                                "in the computation of M0 and M1")

        except (KeyError, AttributeError, TypeError):
            log.warning(f"Unimod entry not found for : {mod}")
    return total_mod_composition


def compute_intensities(df_peptides, unlabelled_aa):
    """Compute isotopologues intensities from peptide sequences.

    Parameters
    ----------
    df_peptides : pandas.Dataframe
        Dataframe with column 'sequence' and 'charge'
    unlabelled_aa : container object
        Container of unlabelled amino acids.

    Returns
    -------
    pandas.Dataframe
        | Dataframe with :
        |                  sequence, mass,
                           formula, formula_X, M0_NC, M1_NC, M0_12C, M1_12C.
    
    Notes
    -----
    | Supports Xtandem's Post-Translational Modification notation (0.4.0).

    """
    log.info("Reading sequences")
    # Remove potential HTML residues from sequences.
    df_peptides["sequence"]= df_peptides["sequence"].str.replace("&gt;", ">", case = False)
    
    # Extract modifications 
    df_peptides["modification"] = df_peptides["sequence"].str.findall(XTANDEM_MOD_PATTERN)
     
    # Remove modifications and capitalize sequence.
    df_peptides["sequence_without_mod"] = df_peptides["sequence"].str.replace(XTANDEM_MOD_PATTERN, "").str.upper()

    # Check that sequences without modifications are real peptide sequences.
    df_peptides["sequence_to_process"], df_peptides["log"] = zip(*df_peptides["sequence_without_mod"].apply(check_amino_acids))

    # Split labelled and unlabelled amino acids from peptide sequences.
    df_peptides["sequence_labelled"], df_peptides["sequence_unlabelled"] = zip(
            *df_peptides["sequence_to_process"].apply(separate_labelled,
                                                      unlabelled_aa=unlabelled_aa))

    log.info("Computing composition and formula")
    # Get composition from modifications
    df_peptides["composition_mod"] = df_peptides["modification"].apply(get_mods_composition)

    # Get Composition from labelled peptide sequence.
    df_peptides["composition_labelled"] = df_peptides["sequence_labelled"].apply(mass.Composition)

    # Get Composition from unlabelled peptide sequence.
    # Unlabelled amino acids are not a real peptide sequence, 
    # hence the 'parsed_sequence' parameter.
    # See https://pyteomics.readthedocs.io/en/latest/api/mass.html
    df_peptides["composition_unlabelled"] = df_peptides["sequence_unlabelled"].apply(lambda x : mass.Composition(parsed_sequence=x))

    # Compute peptide composition
    df_peptides["composition_peptide_neutral"] = (
          df_peptides["composition_labelled"] 
        + df_peptides["composition_unlabelled"]
        + df_peptides["composition_mod"]
    )

    df_peptides["composition_peptide_with_charge"] = (
          df_peptides["composition_peptide_neutral"]
        + df_peptides["charge"].apply(get_charge_composition)
    )

    # Compute peptide composition with X.
    # Carbon atoms from unlabelled peptide are replaced by element X.
    df_peptides["composition_peptide_with_charge_X"] = (
          df_peptides["composition_labelled"]
        + df_peptides["composition_unlabelled"].apply(convert_atom_C_to_X)
        + df_peptides["composition_mod"]
        + df_peptides["charge"].apply(get_charge_composition)
    )

    # Convert formula to string (instead of mass.Composition).
    df_peptides["formula"] = (
        df_peptides["composition_peptide_with_charge"].apply(formula_to_str)
    )
    df_peptides["formula_X"] = (
        df_peptides["composition_peptide_with_charge_X"].apply(formula_to_str)
    )

    # Compute neutral mass
    log.info("Computing neutral mass")
    df_peptides["neutral_mass"] = (
        df_peptides["composition_peptide_neutral"].map(mass.calculate_mass)
    )

    # Add M0 and M1 in normal conditions.
    log.info("Computing M0 and M1")
    # Can use compute_M0_nl with isotopic abundance twice
    df_peptides["M0_NC"] = (
        df_peptides["composition_peptide_with_charge_X"].apply(
                                                        compute_M0_nl,
                                                        a=normal_abundance)
    )
    df_peptides["M1_NC"] = (
        df_peptides["composition_peptide_with_charge_X"].apply(
                                                        compute_M1_nl,
                                                        a=normal_abundance)
    )    
    df_peptides["M0_12C"] = (
        df_peptides["composition_peptide_with_charge_X"].apply(
                                                        compute_M0_nl,
                                                        a=C12_abundance)
    )
    df_peptides["M1_12C"] = (
        df_peptides["composition_peptide_with_charge_X"].apply(
                                                        compute_M1_nl,
                                                        a=C12_abundance)
    )

    return df_peptides.add_prefix('stfi_')


def cli(args=None):
    """Entry point for seq_to_first_iso's CLI.

    Parameters
    ----------
    args : list of str, optional
        CLI arguments, args are used for testing (default is None for CLI).

    Returns
    -------
    None
        Writes a tsv file.

    Notes
    -----
    Main function of the script, for use with CLI.

    """
    if not args:
        args = sys.argv[1:]

    options = user_input(args)

    if options.unlabelled_aa:
        log.info(f"Amino acid with default abundance: {options.unlabelled_aa}")

    log.info("Parsing file")
    df = parse_input_file(option.input_file_name,
                          options.sequence_col_name,
                          options.charge_col_name)
    df = compute_intensities(df, options.unlabelled_aa)

    # Choose output filename.
    if not options.output:
        output_file = option.input_file_name.stem + "_stfi.tsv"
    else:
        output_file = options.output + ".tsv"

    column_of_interest = ["stfi_mass", "stfi_formula", 
                          "stfi_M0_NC", "stfi_M1_NC", 
                          "stfi_M0_12C", "stfi_M1_12C"]
    df["column_of_interest"].to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    cli()  # pragma: no cover
