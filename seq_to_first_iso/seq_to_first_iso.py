"""Main module of seq_to_first_iso.

Provide functions to compute M0 and M1 intensities
with labelled and unlabelled amino acids
for the case of a 99.99 % C[12] enrichment.

The command line interface is also defined here.

Example
-------
Running the script directly
    $ python seq_to_first_iso sequences.tsv sequences charges
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
NATURAL_ABUNDANCE : dict
    Dictionary of isotopic abundances with values taken from MIDAs.
C12_ABUNDANCE : dict
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


USAGE_ERROR = ("Usage: python seq-to-first-iso.py filename.tsv "
               "sequence_column_name charge_column_name "
               "[-o output] [-u aa]"
               )
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

# This variable is obsoleted if an natural element is be named X.
USED_ELEMS = "CHONPSX"

# Columns of interest for export
COLUMNS_OF_INTEREST = ["stfi_neutral_mass",
                       "stfi_formula", "stfi_formula_X",
                       "stfi_M0_NC", "stfi_M1_NC",
                       "stfi_M0_12C", "stfi_M1_12C"]

# Set custom logger.
log = logging.getLogger(__name__)
log_formatter = logging.Formatter("[%(asctime)s] %(levelname)-8s: %(message)s",
                                  "%Y-%m-%d, %H:%M:%S")
log_handler = logging.StreamHandler()
log_handler.setFormatter(log_formatter)
log.addHandler(log_handler)
log.setLevel(logging.INFO)

# Default natural isotopic abundances from MIDAs website:
# https://www.ncbi.nlm.nih.gov/CBBresearch/Yu/midas/index.html .
# X is C with default natural abundance.
NATURAL_ABUNDANCE = {"H[1]": 0.999885, "H[2]": 0.000115,
                     "C[12]": 0.9893, "C[13]": 0.0107,
                     "X[12]": 0.9893, "X[13]": 0.0107,
                     "N[14]": 0.99632, "N[15]": 0.00368,
                     "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                     "S[32]": 0.9493, "S[33]": 0.0076, "S[34]": 0.0429}

C12_ABUNDANCE = dict(NATURAL_ABUNDANCE)
C12_PROPORTION = 0.9999
C12_ABUNDANCE["C[12]"] = C12_PROPORTION
C12_ABUNDANCE["C[13]"] = 1-C12_PROPORTION


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
        description=("Read a tsv file with sequences and charges "
                     "and compute intensity of first isotopologues"
                     ))
    # Input file is required as a positional argument.
    parser.add_argument("input_file_name",
                        type=Path, help="file to parse in .tsv format")
    parser.add_argument("sequence_col_name",
                        type=str, help="column name with sequences")
    parser.add_argument("charge_col_name",
                        type=str, help="column name with charges")

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
                  f"current directory '{options.input_file_name.cwd()}'\n"
                  f"{USAGE_ERROR}"
                  )
        sys.exit()

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


def parse_input_file(filename, sep="\t"):
    r"""Parse input file.

    Parameters
    ----------
    filename : str
        Filename, the file can either just have sequences for each line or
        can have have annotations and sequences with a separator in-between.
    sep : str, optional
        Separator for files with annotations (default is ``\t``).

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    FileNotFoundError
        If the input file is not found.
        Exception chaining is explicitly suppressed (from None).
    Exception
        If the input file cannot be read with pandas.
        Exception chaining is explicitly suppressed (from None).

    """
    if not sep:
        log.warning("Separator is empty, default value '\t' used.")
        sep = "\t"
    try:
        df_input = pd.read_csv(filename, sep=sep)
    except FileNotFoundError:
        log.error(f"File {filename} not found!")
        raise FileNotFoundError(f"File {filename} not found!") from None
    except pd.errors.EmptyDataError:
        log.error(f"Cannot read {filename}!")
        raise Exception(f"Cannot read {filename}!") from None
    log.info(f"Read {filename}")
    return df_input


def filter_input_dataframe(dataframe, sequence_col_name, charge_col_name):
    r"""Filter input file with peptide sequences and charges.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Raw dataframe with all input columns
    sequence_col_name : str
        Name of column with peptide sequences
    charge_col_name : str
        Name of column with peptide charges

    Returns
    -------
    pandas.DataFrame
        | With columns :
        |     - "sequence": peptide sequences.
        |     - "charge": peptide charges.

    Raises
    ------
    KeyError
        If the sequence or charge column is not found.

    """
    line_count, row_count = dataframe.shape
    log.info(f"Found {line_count} lines and {row_count} columns")
    if sequence_col_name not in dataframe.columns:
        log.error(f"Column '{sequence_col_name}' not found in data.")
        raise KeyError(f"Column '{sequence_col_name}' not found in data.")
    if charge_col_name not in dataframe.columns:
        log.error(f"Column '{charge_col_name}' not found in data.")
        raise KeyError(f"Column '{charge_col_name}' not found in data.")
    # Keep only sequences and charges column from original dataframe.
    dataframe = dataframe[[sequence_col_name, charge_col_name]]
    # Rename sequences and charges column to internal naming scheme.
    return dataframe.rename(columns={sequence_col_name: "sequence",
                                     charge_col_name: "charge"})


def check_amino_acids(seq):
    r"""Check elements of a sequence are known amino acids.

    Parameters
    ----------
    seq : str
        Peptide sequence.

    Returns
    -------
    Tuple of two str
        | (sequence, "") if the sequence is composed
        |                of allowed amino acids
        | ("", "Unrecognized amino acids.") if the sequence is composed
        |                                   of unallowed amino acids.

    """
    if not(set(seq) - AMINO_ACIDS) and seq:
        return seq, ""
    return "", "Unrecognized amino acids."


def compute_M0(formula, abundance):
    """Compute intensity of the first isotopologue M0.

    Parameters
    ----------
    formula : pyteomics.mass.Composition
        Chemical formula, as a dict of the number of atoms for each element:
        {element_name: number_of_atoms, ...}.
    abundance : dict
        Dictionary of abundances of isotopes:
        {"element_name[isotope_number]": relative abundance, ..}.

    Returns
    -------
    float
        Value of M0.

    Notes
    -----
    Unused. Use compute_M0_nl instead.

    """
    M0_intensity = (
        abundance["C[12]"]**formula["C"]
        * abundance["H[1]"]**formula["H"]
        * abundance["N[14]"]**formula["N"]
        * abundance["O[16]"]**formula["O"]
        * abundance["S[32]"]**formula["S"]
    )
    return M0_intensity


def compute_M1(formula, abundance):
    """Compute intensity of the second isotopologue M1.

    Parameters
    ----------
    formula : pyteomics.mass.Composition
        Chemical formula, as a dict of the number of atoms for each element:
        {element_name: number_of_atoms, ...}.
    abundance : dict
        Dictionary of abundances of isotopes:
        {"element_name[isotope_number]": relative abundance, ..}.

    Returns
    -------
    float
        Value of M1.

    Notes
    -----
    Unused. Use compute_M1_nl instead.

    """
    M1_intensity = (
        (formula["C"] * abundance["C[12]"]**(formula["C"]-1)
         * abundance["C[13]"]
         * abundance["H[1]"]**formula["H"]
         * abundance["N[14]"]**formula["N"]
         * abundance["O[16]"]**formula["O"]
         * abundance["S[32]"]**formula["S"])

        + (formula["H"] * abundance["C[12]"]**formula["C"]
           * abundance["H[1]"]**(formula["H"]-1) * abundance["H[2]"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**formula["S"])

        + (formula["N"] * abundance["C[12]"]**formula["C"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**(formula["N"]-1) * abundance["N[15]"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**formula["S"])

        + (formula["O"] * abundance["C[12]"]**formula["C"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**(formula["O"]-1) * abundance["O[17]"]
           * abundance["S[32]"]**formula["S"])

        + (formula["S"] * abundance["C[12]"]**formula["C"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**(formula["S"]-1) * abundance["S[33]"])
    )
    return M1_intensity


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


def compute_M0_nl(formula, abundance):
    """Compute intensity of the first isotopologue M0.

    Handle element X with specific abundance.

    Parameters
    ----------
    formula : pyteomics.mass.Composition
        Chemical formula, as a dict of the number of atoms for each element:
        {element_name: number_of_atoms, ...}.
    abundance : dict
        Dictionary of abundances of isotopes:
        {"element_name[isotope_number]": relative abundance, ..}.

    Returns
    -------
    float
        Value of M0.

    Notes
    -----
    X represents C with default isotopic abundance.

    """
    M0_intensity = (
        abundance["C[12]"]**formula["C"]
        * abundance["X[12]"]**formula["X"]
        * abundance["H[1]"]**formula["H"]
        * abundance["N[14]"]**formula["N"]
        * abundance["O[16]"]**formula["O"]
        * abundance["S[32]"]**formula["S"]
    )
    return M0_intensity


def compute_M1_nl(formula, abundance):
    """Compute intensity of the second isotopologue M1.

    Handle element X with specific abundance.

    Parameters
    ----------
    formula : pyteomics.mass.Composition
        Chemical formula, as a dict of the number of atoms for each element:
        {element_name: number_of_atoms, ...}.
    abundance : dict
        Dictionary of abundances of isotopes:
        {"element_name[isotope_number]": relative abundance, ..}.

    Returns
    -------
    float
        Value of M1.

    Notes
    -----
    X represents C with default isotopic abundance.

    """
    M1_intensity = (
        (formula["C"] * abundance["C[12]"]**(formula["C"]-1)
         * abundance["C[13]"]
         * abundance["X[12]"]**formula["X"]
         * abundance["H[1]"]**formula["H"]
         * abundance["N[14]"]**formula["N"]
         * abundance["O[16]"]**formula["O"]
         * abundance["S[32]"]**formula["S"])
        + (formula["X"] * abundance["C[12]"]**formula["C"]
           * abundance["X[12]"]**(formula["X"]-1) * abundance["X[13]"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**formula["S"])
        + (formula["H"] * abundance["C[12]"]**formula["C"]
           * abundance["X[12]"]**formula["X"]
           * abundance["H[1]"]**(formula["H"]-1) * abundance["H[2]"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**formula["S"])
        + (formula["N"] * abundance["C[12]"]**formula["C"]
           * abundance["X[12]"]**formula["X"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**(formula["N"]-1) * abundance["N[15]"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**formula["S"])
        + (formula["O"] * abundance["C[12]"]**formula["C"]
           * abundance["X[12]"]**formula["X"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**(formula["O"]-1) * abundance["O[17]"]
           * abundance["S[32]"]**formula["S"])
        + (formula["S"] * abundance["C[12]"]**formula["C"]
           * abundance["X[12]"]**formula["X"]
           * abundance["H[1]"]**formula["H"]
           * abundance["N[14]"]**formula["N"]
           * abundance["O[16]"]**formula["O"]
           * abundance["S[32]"]**(formula["S"]-1) * abundance["S[33]"])
    )
    return M1_intensity


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


def compute_intensities(df_peptides, unlabelled_aa=[]):
    """Compute isotopologues intensities from peptide sequences.

    Parameters
    ----------
    df_peptides : pandas.DataFrame
        Dataframe with column 'sequence' and 'charge'
    unlabelled_aa : container object
        Container of unlabelled amino acids.

    Returns
    -------
    pandas.DataFrame
        | Dataframe with all computed values, compositions and formulas.

    Notes
    -----
    | Supports Xtandem's Post-Translational Modification notation (0.4.0).

    """
    log.info("Reading sequences.")
    # Remove potential HTML residues from sequences.
    df_peptides["sequence_clean"] = (
        df_peptides["sequence"].str.replace("&gt;", ">", case=False)
    )

    # Extract modifications.
    df_peptides["modification"] = (
        df_peptides["sequence_clean"].str.findall(XTANDEM_MOD_PATTERN)
    )

    # Remove modifications and capitalize sequence.
    df_peptides["sequence_without_mod"] = (
        df_peptides["sequence_clean"]
        .str.replace(XTANDEM_MOD_PATTERN, "")
        .str.upper()
    )

    # Check that sequences without modifications are real peptide sequences.
    df_peptides["sequence_to_process"], df_peptides["log"] = (
        zip(*df_peptides["sequence_without_mod"].apply(check_amino_acids))
    )

    # Split labelled and unlabelled amino acids from peptide sequences.
    df_peptides["sequence_labelled"], df_peptides["sequence_unlabelled"] = zip(
        *df_peptides["sequence_to_process"]
        .apply(separate_labelled, unlabelled_aa=unlabelled_aa)
    )

    log.info("Computing composition and formula.")
    # Get composition from modifications
    df_peptides["composition_mod"] = (
        df_peptides["modification"].apply(get_mods_composition)
    )

    # Get Composition from labelled peptide sequence.
    df_peptides["composition_labelled"] = (
        df_peptides["sequence_labelled"].apply(mass.Composition)
    )

    # Get Composition from unlabelled peptide sequence.
    # Unlabelled amino acids are not a real peptide sequence,
    # hence the 'parsed_sequence' parameter.
    # See https://pyteomics.readthedocs.io/en/latest/api/mass.html
    df_peptides["composition_unlabelled"] = (
        df_peptides["sequence_unlabelled"]
        .apply(lambda x: mass.Composition(parsed_sequence=x))
    )

    # Compute peptide composition.
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
        df_peptides["composition_peptide_with_charge_X"]
        .apply(compute_M0_nl, abundance=NATURAL_ABUNDANCE)
    )
    df_peptides["M1_NC"] = (
        df_peptides["composition_peptide_with_charge_X"]
        .apply(compute_M1_nl, abundance=NATURAL_ABUNDANCE)
    )
    df_peptides["M0_12C"] = (
        df_peptides["composition_peptide_with_charge_X"]
        .apply(compute_M0_nl, abundance=C12_ABUNDANCE)
    )
    df_peptides["M1_12C"] = (
        df_peptides["composition_peptide_with_charge_X"]
        .apply(compute_M1_nl, abundance=C12_ABUNDANCE)
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
        Write a tsv file.

    Notes
    -----
    Main function of the script, for use with CLI.

    """
    if not args:
        args = sys.argv[1:]

    options = user_input(args)
    print(options)
    if options.unlabelled_aa:
        log.info(f"Amino acid with default abundance: {options.unlabelled_aa}")

    log.info("Parsing file")
    df_raw = parse_input_file(options.input_file_name)
    df_filtered = filter_input_dataframe(df_raw,
                                         options.sequence_col_name,
                                         options.charge_col_name)
    df_processed = compute_intensities(df_filtered, options.unlabelled_aa)

    # Choose output filename.
    if not options.output:
        output_file = options.input_file_name.stem + "_stfi.tsv"
    else:
        output_file = options.output + ".tsv"

    # Read original file and append STFI data.
    df_old = pd.read_csv(options.input_file_name, sep="\t")
    df_new = pd.concat([df_old, df_processed[COLUMNS_OF_INTEREST]], axis=1)
    df_new.to_csv(output_file, sep="\t", index=False)


def export_to_knime(df_init, df_processed):
    """Export and merged computed intensities.

    Parameters
    ----------
    df_init : pandas.DataFrame
        Input / initial dataframe.
    df_processed : pandas.DataFrame
        Dataframe with all computed values.

    Returns
    -------
    pandas.DataFrame
        | Dataframe with COLUMNS_OF_INTEREST.

    """
    df_new = pd.concat([df_init, df_processed[COLUMNS_OF_INTEREST]], axis=1)
    return df_new


if __name__ == "__main__":
    cli()  # pragma: no cover.
