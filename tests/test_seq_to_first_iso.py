
"""Test files provided."""

import filecmp
from pathlib import Path

import pandas as pd
from pyteomics import mass
from pyteomics.auxiliary.structures import PyteomicsError
import pytest

import seq_to_first_iso as stfi

data_dir = Path(__file__).parent.resolve().joinpath("test_data")
REL = 1e-6


# Tests for pyteomics.
@pytest.fixture(scope="session")
def get_composition():
    return [("ACDE", {'H': 24, 'C': 15, 'O': 9, 'N': 4, 'S': 1}),
            ("VPKER", {'H': 49, 'C': 27, 'O': 8, 'N': 9}),
            ("ACDEFGHIKLMNPQRSTVWY", {'H': 159, 'C': 107, 'O': 30, 'N': 29, 'S': 2}),
            ("", {}),
            ("ACTG", {'H': 22, 'C': 12, 'O': 6, 'N': 4, 'S': 1}),
            ("LILI", {'H': 46, 'C': 24, 'O': 5, 'N': 4}),
            ("LLLL", {'H': 46, 'C': 24, 'O': 5, 'N': 4})
            ]


def test_composition(get_composition):
    assert mass.Composition("ACDE")
    assert mass.Composition("A") + mass.Composition("C")
    assert mass.Composition(parsed_sequence="ACDE") == {'H': 22, 'C': 15, 'O': 8, 'N': 4, 'S': 1}

    for data in get_composition:
        sequence = data[0]
        expected = data[1]
        assert mass.Composition(sequence) == expected


@pytest.fixture(scope="session")
def get_mass():
    return [("", 0),
            ("ACDE", 436.12639936491996),
            ("VPKER", 627.37040957111),
            ("ACDEFGHIKLMNPQRSTVWY", 2394.1249068251295),
            ("ACTG", 350.1260054421),
            ("LILI", 470.34682059221996),
            ("LLLL", 470.34682059221996),
            ("AAAAAALAAAAAIAAAAAAAAAAAAAGAAAAAAAAAAVIYWWSSEPED", 4108.065244683701),
            ]


def test_calculate_mass(get_mass):
    assert mass.calculate_mass("ACDE") == pytest.approx(436.12639936, REL)
    assert mass.calculate_mass(mass.Composition("ACDE")) == pytest.approx(436.12639936, REL)
    assert mass.calculate_mass(parsed_sequence="ACDE") == pytest.approx(418.115834, REL)
    assert mass.calculate_mass("A") == pytest.approx(89.04767846841, REL)

    for data in get_mass:
        sequence = data[0]
        expected = data[1]
        assert mass.calculate_mass(sequence) == pytest.approx(expected, REL)


# def test_fast_mass(get_mass):
#    assert mass.fast_mass("A")
#    assert mass.fast_mass2("A")
#
#    # fast_mass has a problem with empty sequences.
#    for data in get_mass[1:]:
#        sequence = data[0]
#        expected = data[1]
#        assert mass.fast_mass(sequence) == pytest.approx(expected, REL)
#        assert mass.fast_mass2(sequence) == pytest.approx(expected, REL)


# Tests for seq_to_first_iso.
def test_constants():
    assert stfi.AMINO_ACIDS == set("ACDEFGHIKLMNPQRSTVWY")
    assert stfi.NATURAL_ABUNDANCE == {"H[1]": 0.999885, "H[2]": 0.000115,
                                       "C[12]": 0.9893, "C[13]": 0.0107,
                                       "X[12]": 0.9893, "X[13]": 0.0107,
                                       "N[14]": 0.99632, "N[15]": 0.00368,
                                       "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                                       "S[32]": 0.9493, "S[33]": 0.0076, "S[34]": 0.0429}

    assert stfi.C12_ABUNDANCE["C[12]"] == pytest.approx(0.9999, rel=REL)
    assert stfi.C12_ABUNDANCE["C[13]"] == pytest.approx(0.0001, rel=REL)
    assert {k: v for k, v in stfi.C12_ABUNDANCE.items()
            if k not in ["C[12]", "C[13]"]} == {"H[1]": 0.999885, "H[2]": 0.000115,
                                                "X[12]": 0.9893, "X[13]": 0.0107,
                                                "N[14]": 0.99632, "N[15]": 0.00368,
                                                "O[16]": 0.99757, "O[17]": 0.00038, "O[18]": 0.00205,
                                                "S[32]": 0.9493, "S[33]": 0.0076, "S[34]": 0.0429}


@pytest.fixture(scope="session")
def get_separated_sequence():
    """For C,D and W unlabelled."""
    return [("ABCDEFG", ("ABEFG", "CD")),
            ("ABCDEFGCDCDDC", ("ABEFG", "CDCDCDDC")),
            ("", ("", "")),
            ]


def test_separation(get_separated_sequence):
    test_unlabelled = "CDW"
    for data in get_separated_sequence:
        sequence = data[0]
        expected = data[1]
        assert stfi.separate_labelled(sequence, test_unlabelled)
        result = stfi.separate_labelled(sequence, test_unlabelled)
        assert type(result) is tuple
        assert result == expected


def test_parse_input_file(caplog):
    test_file = data_dir.joinpath("sample_sequence.tsv")
    parser_output = stfi.parse_input_file(test_file, sep="")
    assert "Separator is empty" in caplog.text
    parser_output = stfi.parse_input_file(test_file)
    assert type(parser_output) is pd.DataFrame
    assert parser_output.shape == (8, 3)

    with pytest.raises(FileNotFoundError):
        stfi.parse_input_file("not_a_file")

    empty_file = data_dir.joinpath("sample_empty_file.tsv")
    with pytest.raises(Exception):
        stfi.parse_input_file(empty_file)


def test_filter_input_dataframe(caplog):
    test_file = data_dir.joinpath("sample_sequence.tsv")
    parser_output = stfi.parse_input_file(test_file)
    filter_output = stfi.filter_input_dataframe(parser_output, "sequence", "charge")
    assert type(filter_output) is pd.DataFrame
    assert filter_output.shape == (8, 2)

    bad_file = data_dir.joinpath("bad_sample_sequence.txt")
    with pytest.raises(ValueError):
        parser_output = stfi.parse_input_file(bad_file)
        stfi.parse_input_file(parser_output)


@pytest.fixture(scope="session")
def get_amino_acid_sequences():
    return [("PEPTIDE", ("PEPTIDE", "")),
            ("ACDE", ("ACDE", "")),
            ("VPKER", ("VPKER", "")),
            ("ACDEFGHIKLMNPQRSTVWY", ("ACDEFGHIKLMNPQRSTVWY", "")),
            ("ACTGZ", ("", "Unrecognized amino acids.")),
            ("LILI", ("LILI","")),
            ("AAaag", ("", "Unrecognized amino acids.")),
            ("JOJO", ("", "Unrecognized amino acids.")),
            ]


def test_check_amino_acids(get_amino_acid_sequences):
    for data in get_amino_acid_sequences:
        sequence_input, result = data
        assert stfi.check_amino_acids(sequence_input) == result


def test_deprecated_computation_isotopologue():
    test_composition = mass.Composition("ACDE")
    stfi.compute_M0 = stfi.seq_to_first_iso.compute_M0
    stfi.compute_M1 = stfi.seq_to_first_iso.compute_M1
    assert stfi.compute_M0(test_composition, stfi.NATURAL_ABUNDANCE) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0(test_composition, stfi.C12_ABUNDANCE) == pytest.approx(0.911253268, REL)

    assert stfi.compute_M1(test_composition, stfi.NATURAL_ABUNDANCE) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1(test_composition, stfi.C12_ABUNDANCE) == pytest.approx(0.0277650369575, REL)


def test_computation_isotopologue():
    # Standard formula.
    test_composition = mass.Composition("ACDE")
    assert stfi.compute_M0_nl(test_composition, stfi.NATURAL_ABUNDANCE) == pytest.approx(0.77662382, REL)
    assert stfi.compute_M0_nl(test_composition, stfi.C12_ABUNDANCE) == pytest.approx(0.911253268, REL)
    assert stfi.compute_M1_nl(test_composition, stfi.NATURAL_ABUNDANCE) == pytest.approx(0.1484942353, REL)
    assert stfi.compute_M1_nl(test_composition, stfi.C12_ABUNDANCE) == pytest.approx(0.0277650369575, REL)


def test_string_casting():
    assert stfi.formula_to_str(mass.Composition("ACDE")) == "C15H24O9N4S1"
    assert stfi.formula_to_str(mass.Composition("PEPTIDE")) == "C34H53O15N7"
    assert stfi.formula_to_str(mass.Composition("ACDEFGHIKLMNPQRSTVWY")) == "C107H159O30N29S2"


def test_convert_atom_C_to_X():
    assert stfi.convert_atom_C_to_X("ACDE") == mass.Composition({'H': 24, 'O': 9, 'N': 4, 'S': 1, 'X': 15})
    assert stfi.convert_atom_C_to_X("PEPTIDE") == mass.Composition({'H': 53, 'O': 15, 'N': 7, 'X': 34})
    assert stfi.convert_atom_C_to_X("ACDEFGHIKLMNPQRSTVWY") == mass.Composition({'H': 159, 'O': 30, 'N': 29, 'S': 2, 'X': 107})


@pytest.fixture(scope="session")
def get_charges():
    return [(1, mass.Composition({"H":1})),
            (2, mass.Composition({"H":2})),
            (3, mass.Composition({"H":3})),
            (4, mass.Composition({"H":4})),
            (0, mass.Composition({"H":0})),
            ]

def test_get_charge_composition(get_charges):
    get_charge_composition = stfi.seq_to_first_iso.get_charge_composition
    for charge, target_composition in get_charges:
        assert get_charge_composition(charge) == target_composition


@pytest.fixture(scope="session")
def get_mods():
    return [(["Oxidation"], mass.Composition({"O":1})),
            (["Acetyl", "Phospho"], mass.Composition({'H':3, 'C':2, 'O':4, "P":1})),
            (["Acetyl", "Phospho", "not_mod"], mass.Composition({'H':3, 'C':2, 'O':4, "P":1})),
            ([], mass.Composition()),
            ]


def test_get_mods_composition(get_mods, caplog):
    get_mods_composition = stfi.seq_to_first_iso.get_mods_composition
    get_mods_composition(["not_mod"])
    assert "entry not found" in caplog.text
    # Mod with element not in CHONPSX.
    get_mods_composition(["Heme"])
    assert "(Heme) is not supported" in caplog.text
    for data in get_mods:
        modifications = data[0]
        expected = data[1]
        assert get_mods_composition(modifications) == expected


def test_cli_parser(caplog):
    user_parser = stfi.seq_to_first_iso.user_input
    test_file = data_dir.joinpath("sample_sequence.tsv")

    assert user_parser([str(test_file), "sequence", "charge"])
    with pytest.raises(SystemExit):
        user_parser(["not_a_file", "sequence", "charge"])

    # It seems like a space between -u and A here won't be escaped.
    cli_labelled = user_parser([str(test_file), "sequence", "charge", "-uA,c"])
    assert cli_labelled.unlabelled_aa == ["A", "C"]


def test_main(caplog):
    test_file = data_dir.joinpath("sample_sequence.tsv")

    stfi.seq_to_first_iso.cli([str(test_file), "sequence", "charge", "-uA"])
    assert "Amino acid with default abundance" in caplog.text
    assert Path("sample_sequence_stfi.tsv").is_file()
    assert filecmp.cmp("sample_sequence_stfi.tsv", data_dir.joinpath("sample_sequence_stfi_ref_A.tsv"), shallow=False)
