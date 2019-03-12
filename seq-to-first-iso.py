
"""
Read a file composed of sequences of amino acids on each line and
return sequence, mass, formula, M0 and M1 in normal and C[12] conditions
as a tsv file.

Naming conventions for isotopes follow pyteomics's conventions.
"""

import sys
from pathlib import Path

import pandas as pd
from pyteomics import mass


USAGE_ERROR = "Usage: python seq-to-first-iso.py filename"


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

    # Check if usage is correct.
    if len(sys.argv) < 2:
        exit("Error: not enough arguments  provided\n" + USAGE_ERROR)
    elif not Path.is_file(Path(sys.argv[1])):
        exit("Error: file {} does not exist in current directory {}\n".format(
                sys.argv[1], str(Path.cwd()))
             + USAGE_ERROR)

    filename = sys.argv[1]

    # Obtain a list of sequences as string.
    with open(filename, "r") as filin:
        sequences = []
        # Trim the newline caracters.
        for sequence in filin:
            sequences.append(sequence.strip())

    # Using dataframes to speed up calculations.
    df_peptides = pd.DataFrame({"sequence": sequences})

    # Add mass and formulas of sequences
    df_peptides["mass"] = df_peptides["sequence"].map(mass.calculate_mass)
    df_peptides["formula"] = df_peptides["sequence"].map(mass.Composition)

    # M0 and M1 in normal conditions.
    df_peptides["M0_NC"] = df_peptides["formula"].apply(compute_M0,
                                                        a=isotopic_abundance)
    df_peptides["M1_NC"] = df_peptides["formula"].apply(compute_M1,
                                                        a=isotopic_abundance)

    # M0 and M1 in 99.99 % of C[12].
    prop = 0.9999
    isotopic_abundance["C[12]"] = prop
    isotopic_abundance["C[13]"] = 1-prop
    df_peptides["M0_12C"] = df_peptides["formula"].apply(compute_M0,
                                                         a=isotopic_abundance)
    df_peptides["M1_12C"] = df_peptides["formula"].apply(compute_M1,
                                                         a=isotopic_abundance)

    # Take first 100 rows of df_peptides.
    sample_size = 100
    small_df = df_peptides[:sample_size]

    # Import dataframe to tsv file.
    small_df.to_csv("first-iso.tsv", sep="\t")
